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


class GenomicFeature(NamedTuple):
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
            read.is_unmapped)


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
    parser.add_argument('--input_bed',
                        type=Path,
                        required=True,
                        help='Path to the .bed file with genomic features.')
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
                        help='Min. overlap (in bp) between the read and genomic feature.',
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

    features_df = pd.read_csv(args.input_bed, sep='\t',
                              names=['chromosome', 'start', 'end', 'name', 'score', 'strand'])
    features_by_chrom_and_strand: dict[ChromAndStrand, list[GenomicFeature]] = defaultdict(list)
    for name, chromosome, strand, start, end in zip(features_df['name'],
                                                    features_df['chromosome'],
                                                    features_df['strand'],
                                                    features_df['start'],
                                                    features_df['end']):
        genomic_feature = GenomicFeature(name=name,
                                         chromosome=chromosome,
                                         strand=strand,
                                         start=start,
                                         end=end,
                                         interval=py_interval([start, end]))
        features_by_chrom_and_strand[ChromAndStrand(chromosome=genomic_feature.chromosome,
                                                    strand=genomic_feature.strand)].append(genomic_feature)

    read_counts = {feature.name: 0 for feature in chain.from_iterable(features_by_chrom_and_strand.values())}
    features_by_chrom_and_strand = {key: sorted(value, key=lambda x: x.start) for
                                    key, value in features_by_chrom_and_strand.items()}

    # Check whether features are not overlapping
    for feature_list in features_by_chrom_and_strand.values():
        for feature_1, feature_2 in zip(feature_list, feature_list[1:]):
            if feature_1.end > feature_2.start:
                raise ValueError(
                    f"Overlapping features detected:\n"
                    f"{feature_1.name=}, {feature_1.start=}, {feature_1.end=} \n"
                    f"{feature_2.name=}, {feature_2.start=}, {feature_2.end=} \n"
                    f"Please provide a .bed file with non-overlapping features."
                )

    feature_starts_by_chrom_and_strand = {key: [feature.start for feature in value] for
                                          key, value in features_by_chrom_and_strand.items()}
    feature_ends_by_chrom_and_strand = {key: [feature.end for feature in value] for
                                        key, value in features_by_chrom_and_strand.items()}

    output_unsorted_bam_path = output_folder / 'reads_unsorted.bam'
    output_sorted_bam_path = output_folder / 'reads_sorted.bam'
    output_bed_plus_strand_path = Path(output_folder / 'reads_plus_strand.bed')
    output_bed_minus_strand_path = Path(output_folder / 'reads_minus_strand.bed')

    bam_input = pysam.AlignmentFile(Path(args.input_bam), "rb")

    if create_bam_output:
        bam_output = pysam.AlignmentFile(output_unsorted_bam_path, "wb", template=bam_input)

    if create_bed_output:
        open(output_bed_plus_strand_path, 'w').close()  # Create empty files to append on
        bed_output_plus_strand = open(output_bed_plus_strand_path, 'a')

        open(output_bed_minus_strand_path, 'w').close()
        bed_output_minus_strand = open(output_bed_minus_strand_path, 'a')

    for alignment in generate_alignments(bam_input=bam_input,
                                         paired=args.paired_sequencing,
                                         strandendess_type=strandendess_type,
                                         mapq_threshold=args.mapq_threshold):
        feature_list = features_by_chrom_and_strand.get(ChromAndStrand(alignment.chromosome, alignment.strand))
        if not feature_list:
            continue  # for a case that no features are present on given contig, e.g. introns on MtDNA
        aligned_blocks = py_interval(*(
            alignment.read_1.get_blocks() if alignment.read_2 is None else alignment.read_1.get_blocks() + alignment.read_2.get_blocks()))

        feature_starts = feature_starts_by_chrom_and_strand[ChromAndStrand(alignment.chromosome, alignment.strand)]
        feature_ends = feature_ends_by_chrom_and_strand[ChromAndStrand(alignment.chromosome, alignment.strand)]

        alignment_start = aligned_blocks[0][0]
        alignment_end = aligned_blocks[-1][1]

        feature_index_lower_bound = bisect.bisect_left(feature_ends, alignment_start)
        feature_index_upper_bound = bisect.bisect_right(feature_starts, alignment_end)

        aligned_to_feature = False
        for genomic_feature in feature_list[feature_index_lower_bound:feature_index_upper_bound]:
            interval_intersection = genomic_feature.interval & aligned_blocks
            overlap_length = sum([interval[1] - interval[0] for interval in interval_intersection])

            if overlap_length >= overlap_bp_threshold:
                aligned_to_feature = True
                read_counts[genomic_feature.name] += 1
                if create_bed_output:
                    for read in (alignment.read_1, alignment.read_2):
                        if read is None:
                            continue
                        read_blocks = genomic_feature.interval & py_interval(*read.get_blocks())
                        for block in read_blocks:
                            if alignment.strand == '+':
                                bed_output_plus_strand.write(
                                    f"{alignment.chromosome}\t{int(block[0])}\t{int(block[1])}\n")
                            elif alignment.strand == '-':
                                bed_output_minus_strand.write(
                                    f"{alignment.chromosome}\t{int(block[0])}\t{int(block[1])}\n")

        if aligned_to_feature:
            if create_bam_output:
                bam_output.write(alignment.read_1)
                if alignment.read_2 is not None:
                    bam_output.write(alignment.read_2)

    bam_input.close()

    features_df = features_df.set_index('name', drop=False).drop(columns=['score'])
    features_df['count'] = pd.Series(read_counts)
    features_df.to_csv(output_folder / 'read_counts.tsv', sep='\t', index=False)

    if create_bam_output:
        bam_output.close()
        pysam.sort("-o", str(output_sorted_bam_path), str(output_unsorted_bam_path), catch_stdout=False)
        pysam.index(str(output_sorted_bam_path))
        output_unsorted_bam_path.unlink()
    if create_bed_output:
        bed_output_plus_strand.close()
        bed_output_minus_strand.close()
