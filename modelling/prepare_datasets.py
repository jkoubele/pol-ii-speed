import json
import pickle
from collections import defaultdict
from pathlib import Path
from typing import Optional

import pandas as pd
import torch
from tqdm import tqdm

from pol_ii_model_old import GeneData, GeneDataWithSolution
from read_location_model import estimate_phi


def get_gene_data_with_solution(gene_data: GeneData,
                                design_matrix: torch.tensor,
                                library_size: torch.tensor) -> Optional[GeneDataWithSolution]:
    mask_group_1 = (design_matrix == 0).flatten()
    mask_group_2 = (design_matrix == 1).flatten()

    if gene_data.exon_reads[mask_group_1].sum() == 0 or gene_data.exon_reads[mask_group_2].sum() == 0:
        return None

    mu_1 = float((gene_data.exon_reads[mask_group_1] / library_size[mask_group_1]).mean())
    mu_2 = float((gene_data.exon_reads[mask_group_2] / library_size[mask_group_2]).mean())

    nu_1: dict[str, float] = {}
    nu_2: dict[str, float] = {}

    phi_1: dict[str, float] = {}
    phi_2: dict[str, float] = {}

    introns_to_keep: list[str] = []
    threshold_phi = 0.01
    for intron_name in gene_data.intron_names:
        if gene_data.intron_reads[intron_name][mask_group_1].sum() == 0:
            continue
        if gene_data.intron_reads[intron_name][mask_group_2].sum() == 0:
            continue

        coverage_group_1 = gene_data.coverage_density[intron_name][mask_group_1]

        aggregate_density_1 = (coverage_group_1 * gene_data.intron_reads[intron_name][mask_group_1].unsqueeze(1)).sum(
            dim=0)
        aggregate_density_1 /= aggregate_density_1.sum()

        coverage_group_2 = gene_data.coverage_density[intron_name][mask_group_2]
        aggregate_density_2 = (coverage_group_2 * gene_data.intron_reads[intron_name][mask_group_2].unsqueeze(1)).sum(
            dim=0)
        aggregate_density_2 /= aggregate_density_2.sum()

        intron_phi_1 = estimate_phi(aggregate_density_1.numpy())
        intron_phi_2 = estimate_phi(aggregate_density_2.numpy())

        if not (1 - threshold_phi > intron_phi_1 > threshold_phi):
            continue
        if not (1 - threshold_phi > intron_phi_2 > threshold_phi):
            continue

        nu_1[intron_name] = float(
            (gene_data.intron_reads[intron_name][mask_group_1] / library_size[mask_group_1]).mean())
        nu_2[intron_name] = float(
            (gene_data.intron_reads[intron_name][mask_group_2] / library_size[mask_group_2]).mean())

        phi_1[intron_name] = intron_phi_1
        phi_2[intron_name] = intron_phi_2
        introns_to_keep.append(intron_name)

    if not introns_to_keep:
        return None

    gene_data_filtered = GeneData(exon_reads=gene_data.exon_reads,
                                  intron_names=introns_to_keep,
                                  intron_reads={intron_name: reads for intron_name, reads in
                                                gene_data.intron_reads.items() if intron_name in introns_to_keep},
                                  coverage_density={intron_name: density for intron_name, density in
                                                    gene_data.coverage_density.items() if
                                                    intron_name in introns_to_keep})

    gene_data_with_solution = GeneDataWithSolution(gene_data=gene_data_filtered,
                                                   mu_1=mu_1,
                                                   mu_2=mu_2,
                                                   nu_1=nu_1,
                                                   nu_2=nu_2,
                                                   phi_1=phi_1,
                                                   phi_2=phi_2)
    return gene_data_with_solution


# Prepare toy Drosophila dataset
project_path = Path('/cellfile/datapublic/jkoubele/data_pol_ii/drosophila_mutants')
genome_folder = Path('/cellfile/datapublic/jkoubele/reference_genomes/BDGP6.46')

with open(genome_folder / 'protein_coding_genes.json') as file:
    protein_coding_genes = json.load(file)

exonic_reads_df = pd.read_csv(project_path / 'aggregated_counts' / 'exonic_reads.tsv', sep='\t')
intronic_reads_df = pd.read_csv(project_path / 'aggregated_counts' / 'intron_counts.tsv', sep='\t')
library_size_df = pd.read_csv(project_path / 'aggregated_counts' / 'library_size_factors.tsv', sep='\t')

annotation_df = pd.read_csv(project_path / 'sample_annotation' / 'sample_annotation.tsv', sep='\t')

exonic_reads_df = exonic_reads_df.set_index('gene_id')
intronic_reads_df = intronic_reads_df.set_index('name')
library_size_df = library_size_df.set_index('sample_name')

annotation_df = annotation_df.set_index('sample_name')

assert (annotation_df.index == library_size_df.index).all()
assert (annotation_df.index == exonic_reads_df.columns).all()
assert (annotation_df.index == intronic_reads_df.columns).all()

design_matrix = torch.tensor(annotation_df['genotype'] == 'rp2', dtype=torch.float32).unsqueeze(1)
library_size = torch.tensor(library_size_df['library_size_factor'])
library_size = torch.ones_like(
    library_size)  # TODO: Don't forget to uncomment to actually acount for difference in library size

coverage_df_by_sample: dict[str, pd.DataFrame] = {}
for sample_name in annotation_df.index:
    coverage_df = pd.read_parquet(project_path / 'intron_coverage' / sample_name / 'rescaled_intron_coverages.parquet')
    coverage_df = coverage_df.set_index('intron_name')
    coverage_df_by_sample[sample_name] = coverage_df

gene_to_intron_names = defaultdict(list)
for intron_name in intronic_reads_df.index:
    gene_to_intron_names[intron_name.split('_')[0]].append(intron_name)

data_train: list[GeneData] = []
data_train_with_solutions: list[GeneDataWithSolution] = []

for gene_id, row_exon in tqdm(exonic_reads_df.iterrows(),
                              desc='Preparing gene data',
                              total=len(exonic_reads_df)):
    if row_exon.sum() == 0:
        continue
    intron_names = gene_to_intron_names.get(gene_id)
    if not intron_names:
        continue

    gene_introns_df = intronic_reads_df.loc[intron_names]
    gene_introns_df[gene_introns_df.sum(axis=1) > 0]
    if len(gene_introns_df) == 0:
        continue

    intron_reads: dict[str, torch.tensor] = {}
    coverage_density: dict[str, torch.tensor] = {}
    for intron_name, row_intron in gene_introns_df.iterrows():
        intron_reads[intron_name] = torch.tensor(row_intron, dtype=torch.float32)

        densities = []
        for sample_name in annotation_df.index:
            coverage = coverage_df_by_sample[sample_name].loc[intron_name]
            if coverage.sum() == 0:
                coverage += 1  # put uniform density for introns without any coverage
            densities.append(coverage / coverage.sum())

        coverage_density[intron_name] = torch.tensor(densities, dtype=torch.float32)
    gene_data = GeneData(exon_reads=torch.tensor(row_exon, dtype=torch.float32),
                         intron_names=list(gene_introns_df.index),
                         intron_reads=intron_reads,
                         coverage_density=coverage_density)
    data_train.append(gene_data)

    gene_data_with_solution = get_gene_data_with_solution(gene_data, design_matrix, library_size)
    if gene_data_with_solution is not None:
        data_train_with_solutions.append(gene_data_with_solution)

    if len(data_train_with_solutions) > 100:
        break

output_folder = project_path / 'train_data'
output_folder.mkdir(exist_ok=True)

with open(output_folder / 'data_train_with_solutions.pkl', 'wb') as file:
    pickle.dump({'gene_data_with_solutions': data_train_with_solutions,
                 'design_matrix': design_matrix,
                 'library_size': library_size},
                file)
