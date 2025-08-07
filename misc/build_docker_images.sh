#!/bin/bash

script_directory="$(cd "$(dirname "$0")" && pwd)"
repository_path="$(dirname "$script_directory")"
output_folder="$repository_path/docker_images"

# Create output folder if it doesn't exist
mkdir "$output_folder" -p

# build and save docker images
docker build -t pol_ii_speed_tools "$repository_path/dockerfiles/bioinfo_tools"
docker save -o "$output_folder"/pol_ii_speed_tools.tar pol_ii_speed_tools

docker build -t pol_ii_speed_r "$repository_path/dockerfiles/pol_ii_bioconductor"
docker save -o "$output_folder"/pol_ii_speed_r.tar pol_ii_speed_r

