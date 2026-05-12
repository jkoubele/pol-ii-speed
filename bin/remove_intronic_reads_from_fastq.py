#!/usr/bin/env python3

import argparse
import gzip
from pathlib import Path

import pysam
from Bio import SeqIO
from tqdm import tqdm

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--bam_introns',
                        type=Path,
                        required=True,
                        help='Path to the BAM files with introns.')
    parser.add_argument('--input_fastq',
                        type=Path,
                        required=True,
                        default='Path to the input FASTQ files (compressed as .gz).')
    parser.add_argument('--output_fastq',
                        type=Path,
                        required=True,
                        default='Name (path) of the output FASTQ file.')
    args = parser.parse_args()

    intronic_reads_id: set[str] = set()

    with pysam.AlignmentFile(Path(args.bam_introns), "rb") as bam_input:
        for read in tqdm(bam_input, desc="Reading BAM", mininterval=1):
            intronic_reads_id.add(read.query_name)

    batch_size = 1_000_000
    with gzip.open(Path(args.input_fastq), "rt") as handle_in, open(Path(args.output_fastq), "wt") as handle_out:
        batch = []
        for record in tqdm(SeqIO.parse(handle_in, "fastq"), desc="Filtering reads", mininterval=1):
            if record.id not in intronic_reads_id:
                batch.append(record.format("fastq"))
                if len(batch) >= batch_size:
                    handle_out.write("".join(batch))
                    batch.clear()
        if batch:
            handle_out.write("".join(batch))
