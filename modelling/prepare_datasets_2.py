import json
import pickle
from collections import defaultdict
from pathlib import Path

import pandas as pd
import torch
from tqdm import tqdm

from pol_ii_model import GeneData

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

coverage_df_by_sample: dict[str, pd.DataFrame] = {}
for sample_name in annotation_df.index:
    coverage_df = pd.read_parquet(project_path / 'intron_coverage' / sample_name / 'rescaled_intron_coverages.parquet')
    coverage_df = coverage_df.set_index('intron_name')
    coverage_df_by_sample[sample_name] = coverage_df

gene_to_intron_names = defaultdict(list)
for intron_name in intronic_reads_df.index:
    gene_to_intron_names[intron_name.split('_')[0]].append(intron_name)

data_train: list[GeneData] = []

for gene_id, row_exon in tqdm(exonic_reads_df.iterrows(),
                              desc='Preparing gene data',
                              total=len(exonic_reads_df)):
    if row_exon.sum() == 0:
        continue

    mask_group_1 = (design_matrix == 0).flatten()
    mask_group_2 = (design_matrix == 1).flatten()

    intron_names = gene_to_intron_names.get(gene_id)
    if not intron_names:
        continue

    gene_introns_df = intronic_reads_df.loc[intron_names]
    gene_introns_df = gene_introns_df[gene_introns_df.sum(axis=1) > 0]
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

    if gene_data.exon_reads[mask_group_1].sum() <= 0 or gene_data.exon_reads[mask_group_2].sum() <= 0:
        continue
    introns_to_keep = []
    for intron_name in gene_data.intron_names:
        if gene_data.intron_reads[intron_name][mask_group_1].sum() > 0 and gene_data.intron_reads[intron_name][
            mask_group_2].sum() > 0:
            introns_to_keep.append(intron_name)
    if not introns_to_keep:
        continue
    gene_data_cleaned = GeneData(exon_reads=gene_data.exon_reads,
                                 intron_names=introns_to_keep,
                                 intron_reads={key: value for key, value in gene_data.intron_reads.items() if
                                               key in introns_to_keep},
                                 coverage_density={key: value for key, value in gene_data.coverage_density.items() if
                                                   key in introns_to_keep})

    data_train.append(gene_data_cleaned)

output_folder = project_path / 'train_data'
output_folder.mkdir(exist_ok=True)

with open(output_folder / 'data_train_2_groups.pkl', 'wb') as file:
    pickle.dump({'gene_data': data_train,
                 'design_matrix': design_matrix,
                 'library_size': library_size},
                file)
