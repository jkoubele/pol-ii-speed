#!/bin/bash

#SBATCH --job-name=remove_intronic_reads
#SBATCH --ntasks=10
#SBATCH --mem=20G

# Function to display usage information
usage() {
    echo "Usage: $0 -i <input_folder> -o <output_folder> -b <bam_folder> -s <script_folder>"
    exit 1
}

# Variables to hold arguments
input_folder=""
output_folder=""
bam_folder=""
script_folder=""

# Parse command line arguments
while getopts ":i:o:b:s:" opt; do
    case ${opt} in
        i )
            input_folder=$OPTARG
            ;;
        o )
            output_folder=$OPTARG
            ;;
        b )
            bam_folder=$OPTARG
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
[ -z "$bam_folder" ] || [ -z "$script_folder" ]; then
    echo "Error: Missing mandatory arguments"
    usage
fi

# Create output folder if it doesn't exist
mkdir "$output_folder" -p

# Run docker with script extracting the intronic reads
docker run --rm \
-v "$input_folder":/input_folder \
-v "$output_folder":/output_folder \
-v "$bam_folder":/bam_folder \
-v "$script_folder":/script_folder \
--security-opt seccomp=unconfined \
--init \
bioinfo_tools /bin/sh -c "python3 /script_folder/remove_intronic_reads_from_fastq.py \
--bam_introns /bam_folder/intronic_reads_sorted.bam \
--input_fastq /input_folder/R1.fastq.gz \
--output_fastq /output_folder/R1.fastq; \
pigz -f /output_folder/R1.fastq; \
python3 /script_folder/remove_intronic_reads_from_fastq.py \
--bam_introns /bam_folder/intronic_reads_sorted.bam \
--input_fastq /input_folder/R2.fastq.gz \
--output_fastq /output_folder/R2.fastq; \
pigz -f /output_folder/R2.fastq; \
chmod a+rwX -R /output_folder"
