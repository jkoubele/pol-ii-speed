#!/bin/bash

usage() {
    echo "Usage: $0 [-d <docker_image_path>] [-l <slurm_log_folder>]"
    exit 1
}

script_directory="$(cd "$(dirname "$0")" && pwd)"
repository_path="$(dirname "$script_directory")"

docker_image_path="$repository_path"/docker_images/bioinfo_tools.tar
docker_image_path_2="$repository_path"/docker_images/pol_ii_bioconductor.tar
slurm_log_folder="$repository_path"/slurm_logs

# Parse command line arguments
while getopts ":d:l:" opt; do
    case ${opt} in
        d )
            docker_image_path=$OPTARG
            ;;
        l )
            slurm_log_folder=$OPTARG
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

node_list=$(sinfo -N -h -o "%N" | sort | uniq)

for node in $node_list; do
    echo "Submitting job(s) to reload docker on node: $node"
#    sbatch \
#    --output="$slurm_log_folder"/%x/%j_%x.log \
#    --error="$slurm_log_folder"/%x/%j_%x.err \
#    --nodelist="$node" \
#    "$repository_path"/misc/load_docker_image.sh -d "$docker_image_path"

    sbatch \
    --output="$slurm_log_folder"/%x/%j_%x.log \
    --error="$slurm_log_folder"/%x/%j_%x.err \
    --nodelist="$node" \
    "$repository_path"/misc/load_docker_image.sh -d "$docker_image_path_2"
done