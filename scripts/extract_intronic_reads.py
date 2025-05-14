import argparse
import bisect
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import NamedTuple, Optional

import pandas as pd
import pysam
from interval import interval as py_interval
from tqdm import tqdm


@dataclass
class Alignment:
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

    def to_json_dict(self) -> dict:
        return {'chromosome': self.chromosome,
                'strand': self.strand,
                'start': self.start,
                'end': self.end}


class ChromAndStrand(NamedTuple):
    chromosome: str
    strand: str


def ignore_read(read: pysam.AlignedSegment, mapq_threshold=255) -> bool:
    # Consider filtering by 'N' in read.cigarstring
    return True if ((read.mapping_quality < mapq_threshold) or
                    read.is_secondary or
                    read.is_supplementary or
                    read.is_unmapped or
                    'N' in read.cigarstring) else False


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_bam',
                        default='/cellfile/datapublic/jkoubele/data_pol_ii/celegans_mutants/STAR_output/K002000244_89708//Aligned.sortedByCoord.out.bam')
    parser.add_argument('--intron_bed_file',
                        default='/cellfile/datapublic/jkoubele/reference_genomes/WBcel235/introns_filtered.bed')
    parser.add_argument('--output_folder',
                        default='/cellfile/datapublic/jkoubele/data_pol_ii/celegans_mutants/intronic_reads/test')
    parser.add_argument('--strandendess_type',
                        default='2')
    parser.add_argument('--overlap_bp_threshold',
                        default=10)

    args = parser.parse_args()
    output_folder = Path(args.output_folder)

    introns_df = pd.read_csv(args.intron_bed_file, sep='\t',
                             names=['chromosome', 'start', 'end', 'name', 'score', 'strand'])

    output_folder.mkdir(exist_ok=True, parents=True)

    introns_by_chrom_and_strand: dict[ChromAndStrand, list[Intron]] = defaultdict(list)
    all_introns: list[Intron] = []
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
        all_introns.append(intron)

    read_counts = {intron.name: 0 for intron in all_introns}
    introns_by_chrom_and_strand = {key: sorted(value, key=lambda x: x.start) for
                                   key, value in introns_by_chrom_and_strand.items()}

    # TODO: here we assume that introns are not overlapping.
    # We should explicitely check for that and either raise error, or handle it as a specific case.
    intron_starts_by_chrom_and_strand = {key: [intron.start for intron in value] for
                                         key, value in introns_by_chrom_and_strand.items()}
    intron_ends_by_chrom_and_strand = {key: [intron.end for intron in value] for
                                       key, value in introns_by_chrom_and_strand.items()}

    # input_bam_path = input_folder / 'Aligned.sortedByCoord.out.bam'
    output_unsorted_bam_path = output_folder / 'intronic_reads_unsorted.bam'
    output_sorted_bam_path = output_folder / 'intronic_reads_sorted.bam'
    output_bed_plus_strand_path = Path(output_folder / 'intronic_reads_plus_strand.bed')
    output_bed_minus_strand_path = Path(output_folder / 'intronic_reads_minus_strand.bed')

    bam_input = pysam.AlignmentFile(Path(args.input_bam), "rb")

    # TODO: move to arguments
    paired_sequencing = True
    strandendess_type = args.strandendess_type  # either '1' or '2', eligible for paired sequencing only
    overlap_bp_threshold = int(args.overlap_bp_threshold)
    create_bam_output = True
    create_bed_output = True

    assert strandendess_type in ['1', '2']

    if create_bam_output:
        bam_output = pysam.AlignmentFile(output_unsorted_bam_path, "wb", template=bam_input)

    if create_bed_output:
        open(output_bed_plus_strand_path, 'w').close()  # Create empty files to append on
        bed_output_plus_strand = open(output_bed_plus_strand_path, 'a')

        open(output_bed_minus_strand_path, 'w').close()
        bed_output_minus_strand = open(output_bed_minus_strand_path, 'a')

    if paired_sequencing:
        reads_1: dict[str, pysam.AlignedSegment] = {}
        reads_2: dict[str, pysam.AlignedSegment] = {}

        for read in tqdm(bam_input, mininterval=5):
            if ignore_read(read):
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
                assert False, 'Inconsistent detection of paired reads'

            if read_1.reference_name != read_2.reference_name:
                assert False, 'Paired reads aligned to different chromosomes'

            if strandendess_type == '1':
                strand = '+' if read_1.is_forward else '-'
            elif strandendess_type == '2':
                strand = '+' if read_2.is_forward else '-'
            else:
                assert False, 'Invalid strandendess_type'

            alignment = Alignment(chromosome=read_1.reference_name,
                                  strand=strand,
                                  read_1=read_1,
                                  read_2=read_2)
            # Move to function
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

            if aligned_to_intron:
                if create_bam_output:
                    bam_output.write(alignment.read_1)
                    if alignment.read_2 is not None:
                        bam_output.write(alignment.read_2)
                if create_bed_output:
                    for block in aligned_blocks:
                        if alignment.strand == '+':
                            bed_output_plus_strand.write(f"{alignment.chromosome}\t{int(block[0])}\t{int(block[1])}\n")
                        elif alignment.strand == '-':
                            bed_output_minus_strand.write(f"{alignment.chromosome}\t{int(block[0])}\t{int(block[1])}\n")

    elif not paired_sequencing:
        for read in bam_input:
            alignment = Alignment(chromosome=read.reference_name,
                                  strand='+' if read.is_forward else '-',
                                  read_1=read)

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
