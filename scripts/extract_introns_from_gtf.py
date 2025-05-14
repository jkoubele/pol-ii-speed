import argparse
from pathlib import Path
from collections import defaultdict

import pandas as pd
from pybedtools import BedTool, Interval
from tqdm import tqdm


def extract_genomic_features(gtf_file: Path, output_folder: Path, gtf_source: str, min_length=50) -> None:
    assert gtf_source in ('ensembl', 'gencode'), "gtf_source must be either 'ensembl' or 'gencode'."
    gtf_df = pd.read_csv(gtf_file,
                         header=4,
                         names=['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame',
                                'attribute'],
                         delimiter="\t")

    print(f"Features present in the GTF file: {gtf_df['feature'].unique()}")

    # Processing of the reference genome to obtain introns
    gtf_records = BedTool(gtf_file)
    genes = BedTool(
        [record for record in tqdm(gtf_records, desc='Extracting genes') if record.fields[2] == 'gene']).sort()
    exons = BedTool(
        [record for record in tqdm(gtf_records, desc='Extracting exons') if record.fields[2] == 'exon']).sort()

    if gtf_source == 'ensembl':
        utr_3_prime = BedTool(
            [record for record in tqdm(gtf_records, desc="Extracting three_prime_utr") if
             record.fields[2] == 'three_prime_utr']).sort()
        utr_5_prime = BedTool(
            [record for record in tqdm(gtf_records, desc="Extracting five_prime_utr") if
             record.fields[2] == 'five_prime_utr']).sort()
        utr_3_prime.saveas(output_folder / 'utr_3_prime.bed')
        utr_5_prime.saveas(output_folder / 'utr_5_prime.bed')
        introns = genes.subtract(exons, s=True).subtract(utr_3_prime, s=True).subtract(utr_5_prime, s=True).sort()
    elif gtf_source == 'gencode':
        utr = BedTool(
            [record for record in tqdm(gtf_records, desc="Extracting UTR") if
             record.fields[2] == 'UTR']).sort()
        utr.saveas(output_folder / 'utr.bed')
        introns = genes.subtract(exons, s=True).subtract(utr, s=True).sort()
    else:
        assert False, "Invalid gtf_source."

    introns = BedTool([Interval(chrom=x.chrom, start=x.start, end=x.end,
                                name=x.fields[-1].split()[1][1:-2], score='.',
                                strand=x.strand)
                       for x in tqdm(introns)])
    introns = introns.merge(s=True, c=[4, 5, 6], o='distinct').sort()

    introns_filtered: list[Interval] = []
    gene_occurrences = defaultdict(int)
    for intron in tqdm(introns, desc="Filtering introns"):
        # ignoring introns mapped to more than 1 gene, and too short introns (likely artefacts of annotation)
        if ',' not in intron.name and (intron.end - intron.start) >= min_length:
            gene_occurrences[intron.name] += 1
            intron.name = f"{intron.name}_{gene_occurrences[intron.name]}"
            introns_filtered.append(intron)

    BedTool(introns_filtered).saveas(output_folder / 'introns_filtered.bed')
    genes.saveas(output_folder / 'genes.bed')
    exons.saveas(output_folder / 'exons.bed')
    introns.saveas(output_folder / 'introns.bed')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--output_folder')
    parser.add_argument('--gtf_file', help='Path to the input GTF annotation file.')
    parser.add_argument('--gtf_source', help="Either 'ensembl' or 'gencode'.")
    parser.add_argument('--min_length', default=50, help="Min. intron length")
    args = parser.parse_args()
    extract_genomic_features(gtf_file=Path(args.gtf_file),
                             output_folder=Path(args.output_folder),
                             gtf_source=args.gtf_source,
                             min_length=args.min_length)
