import argparse
import gzip
from pathlib import Path

import pysam
from Bio import SeqIO
from tqdm import tqdm

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--bam_introns',
                        default='/cellfile/datapublic/jkoubele/celegans_mutants/intronic_reads/K002000093_54873/intronic_reads_sorted.bam')
    parser.add_argument('--input_fastq',
                        default='/cellfile/datapublic/jkoubele/celegans_mutants/FASTQ/K002000093_54873/R1.fastq.gz')
    parser.add_argument('--output_fastq',
                        default='/cellfile/datapublic/jkoubele/celegans_mutants/FASTQ_without_intronic_reads/K002000093_54873/R1.fastq')
    args = parser.parse_args()

    intronic_reads_id: set[str] = set()

    with pysam.AlignmentFile(Path(args.bam_introns), "rb") as bam_input:
        for read in tqdm(bam_input, desc="Reading BAM", mininterval=5):
            intronic_reads_id.add(read.query_name)

    with gzip.open(Path(args.input_fastq), "rt") as handle_in, open(Path(args.output_fastq), "wt") as handle_out:
        input_iterator = SeqIO.parse(handle_in, "fastq")
        output_iterator = (record for record in input_iterator if record.id not in intronic_reads_id)
        SeqIO.write(output_iterator, handle_out, "fastq")
