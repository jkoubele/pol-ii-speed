#!/bin/bash

#SBATCH --job-name=salmon_quant
#SBATCH --ntasks=10
#SBATCH --mem=20G

# Function to display usage information
usage() {
    echo "Usage: $0 -i <input_folder> -o <output_folder> -s <salmon_index>"
    exit 1
}

# Variables to hold arguments
input_folder=""
output_folder=""
salmon_index=""

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
            salmon_index=$OPTARG
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
if [ -z "$input_folder" ] || [ -z "$output_folder" ] || [ -z "$salmon_index" ]; then
    echo "Error: Missing mandatory arguments"
    usage
fi

# Create output folder if it doesn't exist
mkdir "$output_folder" -p

# Run Salmon
docker run --rm \
-v "$input_folder":/input_folder \
-v "$output_folder":/output_folder \
-v "$salmon_index":/salmon_index \
--security-opt seccomp=unconfined \
--init \
bioinfo_tools /bin/sh -c "salmon quant \
-i /salmon_index \
-l A \
-1 /input_folder/R1.fastq.gz \
-2 /input_folder/R2.fastq.gz \
-p 8 \
-o \output_folder;
chmod a+rwX -R /output_folder"

