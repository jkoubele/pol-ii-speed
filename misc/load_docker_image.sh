#!/bin/bash
#SBATCH --job-name=load_docker_image
#SBATCH --ntasks=3
#SBATCH --mem=20G

usage() {
    echo "Usage: $0 -d <docker_image_path>"
    exit 1
}

docker_image_path=""

# Parse command line arguments
while getopts ":d:" opt; do
    case ${opt} in
        d )
            docker_image_path=$OPTARG
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
if [ -z "$docker_image_path" ]; then
    echo "Error: Missing mandatory arguments"
    usage
fi

docker load -i "$docker_image_path"
