from pathlib import Path

import numpy as np
import pandas as pd
import torch
from tqdm import tqdm

from pol_ii_speed_modeling.pol_ii_model import GeneData, DatasetMetadata


def load_dataset_metadata(design_matrix_file: Path,
                          library_size_factors_file: Path,
                          lrt_metadata_file: Path,
                          reduced_matrices_folder: Path) -> DatasetMetadata:
    design_matrix_df = pd.read_csv(design_matrix_file, sep='\t')
    library_size_factors_df = pd.read_csv(library_size_factors_file, sep='\t')

    design_matrix_df = design_matrix_df.merge(library_size_factors_df,
                                              left_on='sample',
                                              right_on='sample')
    library_sizes = torch.tensor(design_matrix_df.pop('library_size_factor'),
                                 dtype=torch.float32)
    design_matrix_df = design_matrix_df.set_index('sample')

    lrt_metadata = pd.read_csv(lrt_metadata_file, sep='\t')
    reduced_matrices: dict[str, torch.Tensor] = {}
    for test_id in lrt_metadata['test_id']:
        reduced_matrix_df = pd.read_csv(reduced_matrices_folder / f"{test_id}.tsv", sep='\t').set_index('sample')
        if not all(design_matrix_df.index == reduced_matrix_df.index):
            raise ValueError(
                f"Reduced matrix index {reduced_matrix_df.index} does not equal to design matrix index {design_matrix_df.index}.")
        reduced_matrices[test_id] = torch.tensor(reduced_matrix_df.values, dtype=torch.float32)

    dataset_metadata = DatasetMetadata(design_matrix=torch.tensor(design_matrix_df.values, dtype=torch.float32),
                                       library_sizes=library_sizes,
                                       log_library_sizes=torch.log(library_sizes),
                                       feature_names=design_matrix_df.columns.tolist(),
                                       sample_names=design_matrix_df.index.tolist(),
                                       lrt_metadata=lrt_metadata,
                                       reduced_matrices=reduced_matrices)
    return dataset_metadata


def load_gene_data_list(modeled_genes_file: Path,
                        modeled_introns_file: Path,
                        exon_counts_file: Path,
                        intron_counts_file: Path,
                        isoform_length_factors_file: Path,
                        coverage_folder: Path,
                        sample_names: list[str]) -> list[GeneData]:
    modeled_genes_df = pd.read_csv(modeled_genes_file, sep='\t')
    exon_counts_df = pd.read_csv(exon_counts_file, sep='\t').set_index('gene_id')
    exon_counts_df = exon_counts_df[sample_names]

    intron_counts_df = pd.read_csv(intron_counts_file, sep='\t').set_index('intron_id')
    intron_counts_df = intron_counts_df[sample_names]

    isoform_length_factors_df = pd.read_csv(isoform_length_factors_file, sep='\t').set_index('gene_id')
    isoform_length_factors_df = isoform_length_factors_df[sample_names]

    modeled_introns_df = pd.read_csv(modeled_introns_file, sep='\t')

    gene_to_introns = (
        modeled_introns_df
        .groupby("gene_id")["intron_id"]
        .agg(list)
        .to_dict()
    )

    coverage_df_by_sample: dict[str, pd.DataFrame] = {}
    for sample_name in tqdm(sample_names, desc='Loading coverage data'):
        coverage_df = pd.read_parquet(coverage_folder / f'{sample_name}.parquet')
        coverage_df = coverage_df.set_index('intron_name')
        coverage_df_by_sample[sample_name] = coverage_df

    gene_data_list: list[GeneData] = []

    for gene_id in tqdm(modeled_genes_df['gene_id'], desc='Preparing gene data'):
        exon_row = exon_counts_df.loc[gene_id]

        # Sort introns by position in the gene
        gene_intron_names = sorted(
            gene_to_introns[gene_id],
            key=lambda x: int(x.rsplit('_', maxsplit=1)[1]),
        )

        exon_reads = torch.tensor(exon_row.values, dtype=torch.float32)
        intron_reads = torch.tensor(intron_counts_df.loc[gene_intron_names].values, dtype=torch.float32).T

        coverage = torch.tensor(np.stack([
            np.stack([coverage_df_by_sample[sample].loc[intron].values for intron in gene_intron_names])
            for sample in sample_names
        ]), dtype=torch.float32)

        if not torch.allclose(coverage.sum(axis=2), intron_reads):
            raise ValueError(f'Intron coverage for gene {gene_id} does not sum up to its intron read counts.')

        isoform_length_factors = torch.tensor(isoform_length_factors_df.loc[gene_id].values,
                                              dtype=torch.float32)
        isoform_length_offset = torch.log(isoform_length_factors)

        gene_data = GeneData(gene_name=gene_id,
                             intron_names=gene_intron_names,
                             exon_reads=exon_reads,
                             intron_reads=intron_reads,
                             coverage=coverage,
                             isoform_length_offset=isoform_length_offset)
        gene_data_list.append(gene_data)

    return gene_data_list
