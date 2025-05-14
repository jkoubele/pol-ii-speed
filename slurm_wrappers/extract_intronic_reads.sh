#!/bin/bash

# Function to display usage information
usage() {
    echo "Usage: $0 -i <input_folder> -o <output_folder> -g <genome_folder>"
    exit 1
}

wrapper_script_directory="$(cd "$(dirname "$0")" && pwd)"
repository_path="$(dirname "$wrapper_script_directory")"
slurm_log_folder="$repository_path"/slurm_logs
script_folder="$repository_path"/scripts

# Variables to hold arguments
input_folder=""
output_folder=""
genome_folder=""

# Parse command line arguments
while getopts ":i:o:g:" opt; do
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
if [ -z "$input_folder" ] || [ -z "$output_folder" ] || [ -z "$genome_folder" ]; then
    echo "Error: Missing mandatory arguments"
    usage
fi

# Create output folder if it doesn't exist
mkdir "$output_folder" -p

for sub_folder in "$input_folder"/*; do
  sample_name=$(basename "$sub_folder")
  echo "Submitting sample $sample_name"
  sbatch \
  --output="$slurm_log_folder"/%x/%j_%x.log \
  --error="$slurm_log_folder"/%x/%j_%x.err \
  "$repository_path"/jobs/extract_intronic_reads.sh \
  -i "$sub_folder" \
  -o "$output_folder"/"$sample_name" \
  -g "$genome_folder" \
  -s "$script_folder"
done
