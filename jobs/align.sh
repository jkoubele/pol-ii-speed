#!/bin/bash

#SBATCH --job-name=align
#SBATCH --ntasks=15
#SBATCH --mem=80G

# Function to display usage information
usage() {
    echo "Usage: $0 -i <input_folder> -o <output_folder> -s <star_index>"
    exit 1
}

# Variables to hold arguments
input_folder=""
output_folder=""
star_index=""

# Parse command line arguments
while getopts ":i:o:s:" opt; do
    case ${opt} in
        i )
            input_folder=$OPTARG
            ;;
        o )
            output_folder=$OPTARG
            ;;
        s )
            star_index=$OPTARG
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
if [ -z "$input_folder" ] || [ -z "$output_folder" ] || [ -z "$star_index" ]; then
    echo "Error: Missing mandatory arguments"
    usage
fi

# Create output folder if it doesn't exist
mkdir "$output_folder" -p

# Run STAR aligner
docker run --rm \
-v "$input_folder":/input_folder \
-v "$output_folder":/output_folder \
-v "$star_index":/star_index \
--security-opt seccomp=unconfined \
--init \
bioinfo_tools /bin/sh -c "STAR \
--runThreadN 13 \
--genomeDir /star_index \
--readFilesIn /input_folder/R1.fastq.gz /input_folder/R2.fastq.gz \
--readFilesCommand zcat \
--outFileNamePrefix /output_folder/ \
--outSAMtype BAM SortedByCoordinate \
--outSAMattributes All \
--quantMode GeneCounts \
--limitBAMsortRAM 50000000000 \
--peOverlapNbasesMin 10; \
samtools index /output_folder/Aligned.sortedByCoord.out.bam; \
chmod a+rwX -R /output_folder"
