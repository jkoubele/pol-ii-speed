#!/bin/bash

#SBATCH --job-name=compute_coverage
#SBATCH --ntasks=10
#SBATCH --mem=20G

# Function to display usage information
usage() {
    echo "Usage: $0 -i <input_folder> -o <output_folder> -g <genome_folder> -f <fai_file_name> -s <script_folder>"
    exit 1
}

# Variables to hold arguments
input_folder=""
output_folder=""
genome_folder=""
fai_file_name=""
script_folder=""

# Parse command line arguments
while getopts ":i:o:g:f:s:" opt; do
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
        f )
            fai_file_name=$OPTARG
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
if [ -z "$input_folder" ] || [ -z "$output_folder" ] || \
[ -z "$genome_folder" ] || [ -z "$fai_file_name" ] || [ -z "$script_folder" ]; then
    echo "Error: Missing mandatory arguments"
    usage
fi

# Create output folder if it doesn't exist
mkdir "$output_folder" -p

# Compute coverage
docker run --rm \
-v "$input_folder":/input_folder \
-v "$output_folder":/output_folder \
-v "$genome_folder":/genome_folder \
-v "$script_folder":/script_folder \
--security-opt seccomp=unconfined \
--init \
pol_ii_bioconductor /bin/sh -c "
for orientation in plus minus; do
  bedtools genomecov -bga -split \
  -i /input_folder/intronic_reads_\${orientation}_strand.bed.gz \
  -g /genome_folder/$fai_file_name \
  > /output_folder/coverage_\${orientation}_strand.bedGraph; \
  pigz -f /output_folder/coverage_\${orientation}_strand.bedGraph; \
done; \
Rscript /script_folder/compute_rescaled_coverage.R \
--bed_graph_plus /output_folder/coverage_plus_strand.bedGraph.gz \
--bed_graph_minus /output_folder/coverage_minus_strand.bedGraph.gz \
--introns_file /genome_folder/introns_filtered.bed \
--output_folder /output_folder; \
chmod a+rwX -R /output_folder"

