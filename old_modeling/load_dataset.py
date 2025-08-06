from collections import defaultdict
from enum import StrEnum
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from tqdm import tqdm

from modeling import GeneData, DatasetMetadata

MIN_SAMPLES_WITH_NONZERO_READ_COUNT = 3


class IgnoringGeneReasons(StrEnum):
    NOT_IN_EXON_READS = 'Gene not found in the table with exonic reads.'
    NOT_ENOUGH_EXON_READS = f'Gene has less than {MIN_SAMPLES_WITH_NONZERO_READ_COUNT} samples with non-zero read counts.'
    NO_INTRONS_IN_TABLE = 'Gene has no introns in the intronic count matrix.'
    NO_INTRONS_AFTER_FILTERING = 'Gene has no introns that would have enough non-zero read counts.'


class IgnoringIntronReasons(StrEnum):
    NOT_ENOUGH_INTRON_READS = f'Intron has less than {MIN_SAMPLES_WITH_NONZERO_READ_COUNT} samples with non-zero read count.'


def load_dataset_matedata(design_matrix_file: Path,
                          library_size_factors_file: Path) -> DatasetMetadata:
    design_matrix_df = pd.read_csv(design_matrix_file)
    library_size_factors_df = pd.read_csv(library_size_factors_file, sep='\t')

    design_matrix_df = design_matrix_df.merge(library_size_factors_df,
                                              left_on='sample',
                                              right_on='sample')
    library_sizes = torch.tensor(design_matrix_df.pop('library_size_factor'),
                                 dtype=torch.float32)
    design_matrix_df = design_matrix_df.set_index('sample')
    
    dataset_metadata = DatasetMetadata(design_matrix=torch.tensor(design_matrix_df.values, dtype=torch.float32),
                                       library_sizes=library_sizes,
                                       log_library_sizes=torch.log(library_sizes),
                                       feature_names=design_matrix_df.columns.tolist(),
                                       sample_names=design_matrix_df.index.tolist())
    return dataset_metadata

 
def load_gene_data_list(gene_names_file: Path,
                        exon_counts_file: Path,
                        intron_counts_file: Path,
                        coverage_folder: Path,
                        sample_names: list[str],
                        log_output_folder: Path) -> list[GeneData]:
    log_output_folder.mkdir(exist_ok=True, parents=True)
    gene_names_df = pd.read_csv(gene_names_file)

    exon_counts_df = pd.read_csv(exon_counts_file, sep='\t').set_index('gene_id')
    exon_counts_df = exon_counts_df[sample_names]

    intron_counts_df = pd.read_csv(intron_counts_file, sep='\t').set_index('intron_id')
    intron_counts_df = intron_counts_df[sample_names]

    gene_to_intron_names = defaultdict(list)
    for intron_name in intron_counts_df.index:
        gene_to_intron_names[intron_name.rsplit('_', maxsplit=1)[0]].append(intron_name)

    coverage_df_by_sample: dict[str, pd.DataFrame] = {}
    for sample_name in tqdm(sample_names, desc='Loading coverage data'):
        coverage_df = pd.read_parquet(coverage_folder / f'{sample_name}.parquet')
        coverage_df = coverage_df.set_index('intron_name')
        coverage_df_by_sample[sample_name] = coverage_df

    # Check for missing samples in count matrices    
    missing_exon_samples = set(sample_names) - set(exon_counts_df.columns)
    if missing_exon_samples:
        raise ValueError(f"The following samples are missing in exon_counts_df: {sorted(missing_exon_samples)}")
    missing_intron_samples = set(sample_names) - set(intron_counts_df.columns)
    if missing_intron_samples:
        raise ValueError(f"The following samples are missing in intron_counts_df: {sorted(missing_intron_samples)}")

    gene_data_list: list[GeneData] = []

    ignored_genes: dict[str, str] = {}
    ignored_introns: dict[str, str] = {}

    for gene_id in tqdm(gene_names_df['gene_id'].unique(), desc='Preparing gene data'):
        if gene_id not in exon_counts_df.index:
            ignored_genes[gene_id] = IgnoringGeneReasons.NOT_IN_EXON_READS.value
            continue

        exon_row = exon_counts_df.loc[gene_id]
        if (exon_row > 0).sum() < MIN_SAMPLES_WITH_NONZERO_READ_COUNT:
            ignored_genes[gene_id] = IgnoringGeneReasons.NOT_ENOUGH_EXON_READS.value
            continue

        gene_intron_names = gene_to_intron_names[gene_id]
        if not gene_intron_names:
            ignored_genes[gene_id] = IgnoringGeneReasons.NO_INTRONS_IN_TABLE.value
            continue

        introns_to_keep: list[str] = []
        for intron_name in gene_intron_names:
            if (intron_counts_df.loc[intron_name] > 0).sum() < MIN_SAMPLES_WITH_NONZERO_READ_COUNT:
                ignored_introns[intron_name] = IgnoringIntronReasons.NOT_ENOUGH_INTRON_READS.value
                continue
            else:
                introns_to_keep.append(intron_name)

        if not introns_to_keep:
            ignored_genes[gene_id] = IgnoringGeneReasons.NO_INTRONS_AFTER_FILTERING.value
            continue

        # Sort by intron position in the gene
        introns_to_keep = sorted(introns_to_keep, key=lambda x: int(x.rsplit('_', maxsplit=1)[1]))

        exon_reads = torch.tensor(exon_row.values, dtype=torch.float32)
        intron_reads = torch.tensor(intron_counts_df.loc[introns_to_keep].values, dtype=torch.float32).T

        coverage_of_bases = torch.tensor(np.stack([
            np.stack([coverage_df_by_sample[sample].loc[intron].values for intron in introns_to_keep])
            for sample in sample_names
        ]), dtype=torch.float32)

        coverage_sum = coverage_of_bases.sum(axis=2).unsqueeze(2)

        coverage_density = torch.where(
            coverage_sum > 0,
            coverage_of_bases / coverage_sum,
            torch.zeros_like(coverage_of_bases)
        )

        read_coverage = coverage_density * intron_reads.unsqueeze(2)

        gene_data = GeneData(gene_name=gene_id,
                             intron_names=introns_to_keep,
                             exon_reads=exon_reads,
                             intron_reads=intron_reads,
                             coverage=read_coverage)
        gene_data_list.append(gene_data)

    ignored_genes_df = pd.DataFrame(data={'gene': list(ignored_genes.keys()),
                                          'reason': list(ignored_genes.values())})
    ignored_introns_df = pd.DataFrame(data={'intron': list(ignored_introns.keys()),
                                            'reason': list(ignored_introns.values())})
    ignored_genes_df.to_csv(log_output_folder / 'ignored_genes.csv', index=False)
    ignored_introns_df.to_csv(log_output_folder / 'ignored_introns.csv', index=False)
    
    return gene_data_list
    

def load_dataset_from_results_folder(results_folder: Path,
                                     log_output_folder: Path) -> tuple[list[GeneData], DatasetMetadata]:    
    design_matrix_file = results_folder / 'design_matrix' / 'design_matrix.csv'
    library_size_factors_file = results_folder / 'aggregated_counts' / 'library_size_factors.tsv' 
    
    dataset_metadata = load_dataset_matedata(design_matrix_file, library_size_factors_file)
    
    # Another function arguments
    log_output_folder = output_folder / 'logs'
    gene_names_file = results_folder / 'gene_names' / 'protein_coding_genes.csv'
    exon_counts_file = results_folder / 'aggregated_counts' / 'exon_counts.tsv'
    intron_counts_file = results_folder / 'aggregated_counts' / 'intron_counts.tsv'
    coverage_folder = results_folder / 'rescaled_coverage'
    
    gene_data_list = load_gene_data_list(gene_names_file=gene_names_file,
                        exon_counts_file=exon_counts_file, 
                        intron_counts_file=intron_counts_file, 
                        coverage_folder=coverage_folder,
                        sample_names = dataset_metadata.sample_names,
                        log_output_folder=log_output_folder                        
                        )
    
    return dataset_metadata, gene_data_list
    
    
    


if __name__ == "__main__":    
    results_folder = Path('/cellfile/datapublic/jkoubele/drosophila_mutants/results')
    output_folder =  Path('/cellfile/datapublic/jkoubele/drosophila_mutants/results/chunk_model_results/test_chunk')
    
    dataset_metadata, gene_data_list = load_dataset_from_results_folder(results_folder, output_folder)    
   
