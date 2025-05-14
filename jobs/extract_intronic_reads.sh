#!/bin/bash

#SBATCH --job-name=extract_intronic_reads
#SBATCH --ntasks=10
#SBATCH --mem=25G

# Function to display usage information
usage() {
    echo "Usage: $0 -i <input_folder> -o <output_folder> -g <genome_folder> -s <script_folder>"
    exit 1
}

# Variables to hold arguments
input_folder=""
output_folder=""
genome_folder=""
script_folder=""

# Parse command line arguments
while getopts ":i:o:g:s:" opt; do
    case ${opt} in
        i )
            input_folder=$OPTARG
            ;;
        o )
            output_folder=$OPTARG
            ;;
        g )
            genome_folder=$OPTARG
            ;;
        s )
            script_folder=$OPTARG
            ;;
        \? )
            echo "Invalid option: $OPTARG" 1>&2
            usage
            ;;
        : )
            echo "Invalid option: $OPTARG requires an argument" 1>&2
            usage
            ;;
    esac
done

# Check if mandatory arguments are provided
if [ -z "$input_folder" ] || [ -z "$output_folder" ] ||
[ -z "$genome_folder" ] || [ -z "$script_folder" ]; then
    echo "Error: Missing mandatory arguments"
    usage
fi

# Create output folder if it doesn't exist
mkdir "$output_folder" -p

# Extract intronic reads
docker run --rm \
-v "$input_folder":/input_folder \
-v "$output_folder":/output_folder \
-v "$genome_folder":/genome_folder \
-v "$script_folder":/script_folder \
--security-opt seccomp=unconfined \
--init \
bioinfo_tools /bin/sh -c "python3 /script_folder/extract_intronic_reads.py \
--input_bam /input_folder/Aligned.sortedByCoord.out.bam \
--intron_bed_file /genome_folder/introns_filtered.bed \
--output_folder /output_folder; \
sort -k 1,1 -k 2,2n /output_folder/intronic_reads_plus_strand.bed > /output_folder/tmp_plus_strand.bed; \
mv /output_folder/tmp_plus_strand.bed /output_folder/intronic_reads_plus_strand.bed; \
sort -k 1,1 -k 2,2n /output_folder/intronic_reads_minus_strand.bed > /output_folder/tmp_minus_strand.bed; \
mv /output_folder/tmp_minus_strand.bed /output_folder/intronic_reads_minus_strand.bed; \
pigz -f /output_folder/intronic_reads_plus_strand.bed; \
pigz -f /output_folder/intronic_reads_minus_strand.bed; \
chmod a+rwX -R /output_folder"
