#!/bin/bash

script_directory="$(cd "$(dirname "$0")" && pwd)"
repository_path="$(dirname "$script_directory")"
output_folder="$repository_path/docker_images"

# Create output folder if it doesn't exist
mkdir "$output_folder" -p

# Ensure that latest debian base image is used
docker pull debian:latest

# build and save docker images
docker build -t bioinfo_tools "$repository_path/dockerfiles/bioinfo_tools"
docker save -o "$output_folder"/bioinfo_tools.tar bioinfo_tools
