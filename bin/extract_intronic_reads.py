#!/usr/bin/env python3

import argparse
import bisect
from collections import defaultdict
from enum import StrEnum
from itertools import chain
from pathlib import Path
from typing import NamedTuple, Optional, Iterator

import pandas as pd
import pysam
from interval import interval as py_interval
from tqdm import tqdm


class StrandednessType(StrEnum):
    FORWARD = 'forward'
    REVERSE = 'reverse'


class Alignment(NamedTuple):
    chromosome: str
    strand: str
    read_1: pysam.AlignedSegment
    read_2: Optional[pysam.AlignedSegment] = None


class Intron(NamedTuple):
    name: str
    chromosome: str
    strand: str
    start: int
    end: int
    interval: py_interval


class ChromAndStrand(NamedTuple):
    chromosome: str
    strand: str


def ignore_read(read: pysam.AlignedSegment, mapq_threshold=255) -> bool:
    return (read.mapping_quality < mapq_threshold or
            read.is_secondary or
            read.is_supplementary or
            read.is_unmapped or
            'N' in read.cigarstring)


def generate_alignments(bam_input: pysam.AlignmentFile,
                        paired: bool,
                        strandendess_type: StrandednessType,
                        mapq_threshold: int = 255) -> Iterator[Alignment]:
    if not paired:
        raise NotImplementedError("Unpaired mode not implemented yet :( ")

    reads_1: dict[str, pysam.AlignedSegment] = {}
    reads_2: dict[str, pysam.AlignedSegment] = {}

    for read in tqdm(bam_input, mininterval=1):
        if ignore_read(read, mapq_threshold=mapq_threshold):
            continue

        if read.is_read1 and read.query_name not in reads_2:
            reads_1[read.query_name] = read
            continue
        elif read.is_read2 and read.query_name not in reads_1:
            reads_2[read.query_name] = read
            continue

        if read.is_read1 and read.query_name in reads_2:
            read_1 = read
            read_2 = reads_2.pop(read.query_name)
        elif read.is_read2 and read.query_name in reads_1:
            read_1 = reads_1.pop(read.query_name)
            read_2 = read
        else:
            raise ValueError('Inconsistent detection of paired reads!')

        if read_1.reference_name != read_2.reference_name:
            raise ValueError('Paired reads aligned to different chromosomes!')

        if strandendess_type == StrandednessType.FORWARD:
            strand = '+' if read_1.is_forward else '-'
        elif strandendess_type == StrandednessType.REVERSE:
            strand = '+' if read_2.is_forward else '-'
        else:
            raise RuntimeError(f'Invalid strandendess {strandendess_type}')

        yield Alignment(chromosome=read_1.reference_name,
                        strand=strand,
                        read_1=read_1,
                        read_2=read_2)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_bam',
                        type=Path,
                        required=True,
                        help='Path to the input .bam file.')
    parser.add_argument('--intron_bed_file',
                        type=Path,
                        required=True,
                        help='Path to the .bed file with introns.')
    parser.add_argument('--output_folder',
                        type=Path,
                        default='.')
    parser.add_argument('--strandedness',
                        type=StrandednessType,
                        choices=[StrandednessType.FORWARD, StrandednessType.REVERSE],
                        default=StrandednessType.REVERSE,
                        help="Either 'forward' or 'reverse'."
                        )
    parser.add_argument('--overlap_bp_threshold',
                        type=int,
                        help='Min. overlap with intron needed to classify the read as intronic.',
                        default=10)

    parser.add_argument('--mapq_threshold',
                        type=int,
                        help='Min. MAPQ value used to filter reads in the .bam file. '
                             'The default value of 255 corresponds to uniquely mapped reads by STAR. '
                             'Please note that different alignment tools may use different values of MAPQ '
                             'for uniquely mapped reads, e.g. HISAT2 assigns MAPQ of 60 to them.',
                        default=255)

    parser.add_argument('--unpaired_sequencing', action='store_false', dest='paired_sequencing',
                        help='Use for unpaired sequencing (paired sequencing is assumed otherwise).')

    parser.add_argument('--no_bam_output', action='store_false', dest='create_bam_output',
                        help='Disable .bam output (enabled by default).')

    parser.add_argument('--no_bed_output', action='store_false', dest='create_bed_output',
                        help='Disable .bed output (enabled by default).')

    args = parser.parse_args()

    output_folder = args.output_folder
    mapq_threshold = args.mapq_threshold
    paired_sequencing = args.paired_sequencing
    strandendess_type = args.strandedness
    overlap_bp_threshold = args.overlap_bp_threshold
    create_bam_output = args.create_bam_output
    create_bed_output = args.create_bed_output

    output_folder.mkdir(exist_ok=True, parents=True)

    introns_df = pd.read_csv(args.intron_bed_file, sep='\t',
                             names=['chromosome', 'start', 'end', 'name', 'score', 'strand'])
    introns_by_chrom_and_strand: dict[ChromAndStrand, list[Intron]] = defaultdict(list)
    for name, chromosome, strand, start, end in zip(introns_df['name'],
                                                    introns_df['chromosome'],
                                                    introns_df['strand'],
                                                    introns_df['start'],
                                                    introns_df['end']):
        intron = Intron(name=name,
                        chromosome=chromosome,
                        strand=strand,
                        start=start,
                        end=end,
                        interval=py_interval([start, end]))
        introns_by_chrom_and_strand[ChromAndStrand(chromosome=intron.chromosome,
                                                   strand=intron.strand)].append(intron)

    read_counts = {intron.name: 0 for intron in chain.from_iterable(introns_by_chrom_and_strand.values())}
    introns_by_chrom_and_strand = {key: sorted(value, key=lambda x: x.start) for
                                   key, value in introns_by_chrom_and_strand.items()}

    # Check whether introns are not overlapping
    for intron_list in introns_by_chrom_and_strand.values():
        for intron_1, intron_2 in zip(intron_list, intron_list[1:]):
            if intron_1.end > intron_2.start:
                raise ValueError(
                    f"Overlapping introns detected:\n"
                    f"{intron_1.name=}, {intron_1.start=}, {intron_1.end=} \n"
                    f"{intron_2.name=}, {intron_2.start=}, {intron_2.end=} \n"
                    f"Please provide a .bed file with non-overlapping introns."
                )

    intron_starts_by_chrom_and_strand = {key: [intron.start for intron in value] for
                                         key, value in introns_by_chrom_and_strand.items()}
    intron_ends_by_chrom_and_strand = {key: [intron.end for intron in value] for
                                       key, value in introns_by_chrom_and_strand.items()}

    output_unsorted_bam_path = output_folder / 'intronic_reads_unsorted.bam'
    output_sorted_bam_path = output_folder / 'intronic_reads_sorted.bam'
    output_bed_plus_strand_path = Path(output_folder / 'intronic_reads_plus_strand.bed')
    output_bed_minus_strand_path = Path(output_folder / 'intronic_reads_minus_strand.bed')

    bam_input = pysam.AlignmentFile(Path(args.input_bam), "rb")

    if create_bam_output:
        bam_output = pysam.AlignmentFile(output_unsorted_bam_path, "wb", template=bam_input)

    if create_bed_output:
        open(output_bed_plus_strand_path, 'w').close()  # Create empty files to append on
        bed_output_plus_strand = open(output_bed_plus_strand_path, 'a')

        open(output_bed_minus_strand_path, 'w').close()
        bed_output_minus_strand = open(output_bed_minus_strand_path, 'a')

    for alignment in generate_alignments(bam_input=bam_input,
                                         paired=True,
                                         strandendess_type=strandendess_type,
                                         mapq_threshold=args.mapq_threshold):
        intron_list = introns_by_chrom_and_strand.get(ChromAndStrand(alignment.chromosome, alignment.strand))
        if not intron_list:
            continue  # for a case that no introns are present on given contig, e.g. for MtDNA
        aligned_blocks = py_interval(*(
            alignment.read_1.get_blocks() if alignment.read_2 is None else alignment.read_1.get_blocks() + alignment.read_2.get_blocks()))

        intron_starts = intron_starts_by_chrom_and_strand[ChromAndStrand(alignment.chromosome, alignment.strand)]
        intron_ends = intron_ends_by_chrom_and_strand[ChromAndStrand(alignment.chromosome, alignment.strand)]

        alignment_start = aligned_blocks[0][0]
        alignment_end = aligned_blocks[-1][1]

        intron_index_lower_bound = bisect.bisect_left(intron_ends, alignment_start)
        intron_index_upper_bound = bisect.bisect_right(intron_starts, alignment_end)

        aligned_to_intron = False
        for intron in intron_list[intron_index_lower_bound:intron_index_upper_bound]:
            interval_intersection = intron.interval & aligned_blocks
            overlap_length = sum([interval[1] - interval[0] for interval in interval_intersection])

            if overlap_length >= overlap_bp_threshold:
                aligned_to_intron = True
                read_counts[intron.name] += 1
                if create_bed_output:
                    for read in (alignment.read_1, alignment.read_2):
                        if read is None:
                            continue
                        read_blocks = intron.interval & py_interval(*read.get_blocks())
                        for block in read_blocks:
                            if alignment.strand == '+':
                                bed_output_plus_strand.write(
                                    f"{alignment.chromosome}\t{int(block[0])}\t{int(block[1])}\n")
                            elif alignment.strand == '-':
                                bed_output_minus_strand.write(
                                    f"{alignment.chromosome}\t{int(block[0])}\t{int(block[1])}\n")

        if aligned_to_intron:
            if create_bam_output:
                bam_output.write(alignment.read_1)
                if alignment.read_2 is not None:
                    bam_output.write(alignment.read_2)

    bam_input.close()

    introns_df = introns_df.set_index('name', drop=False).drop(columns=['score'])
    introns_df['count'] = pd.Series(read_counts)
    introns_df.to_csv(output_folder / 'intron_read_counts.tsv', sep='\t', index=False)

    if create_bam_output:
        bam_output.close()
        pysam.sort("-o", str(output_sorted_bam_path), str(output_unsorted_bam_path), catch_stdout=False)
        pysam.index(str(output_sorted_bam_path))
        output_unsorted_bam_path.unlink()
    if create_bed_output:
        bed_output_plus_strand.close()
        bed_output_minus_strand.close()
