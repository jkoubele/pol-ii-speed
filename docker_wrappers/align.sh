#!/bin/bash

#SBATCH --job-name=align
#SBATCH --ntasks=15
#SBATCH --mem=60G

# Function to display usage information
usage() {
    echo "Usage: $0 -i <input_folder> -o <output_folder> -g <genome_folder>"
    exit 1
}

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

# Run STAR aligner
docker run --rm -v "$input_folder":/input_folder -v "$output_folder":/output_folder \
-v "$genome_folder":/genome_folder --security-opt seccomp=unconfined \
bioinfo_tools /bin/sh -c "STAR \
--runThreadN 14 \
--genomeDir /genome_folder/STAR_index \
--readFilesIn /input_folder/R1.fastq.gz /input_folder/R2.fastq.gz \
--readFilesCommand zcat \
--outFileNamePrefix /output_folder/ \
--outSAMtype BAM SortedByCoordinate \
--outSAMattributes GX GN \
--quantMode GeneCounts \
--limitBAMsortRAM 50000000000 \
--peOverlapNbasesMin 10; \
samtools index /output_folder/Aligned.sortedByCoord.out.bam; \
chmod 777 -R /output_folder"
