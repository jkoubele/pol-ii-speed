# pol-ii-speed
Repo with a pipeline using total RNA-seq to estimate the changes of RNA Polymerase II
elongation speed, and also the changes in the speed of intron splicing.

The code is organized as a Nextflow pipeline. To run it, you will need FASTQ
data of total RNA-seq; please read the documentation below for more details.

The model itself is implemented using PyTorch; the code can be found in the [pol_ii_speed_modeling](pol_ii_speed_modeling) folder, 
which can be also installed as a pip package.

## Running the pipeline
Besides the FASTQ data, you will need to prepare two files to run the pipeline:
 * ```samplesheet.csv```, containing sample annotation.
 * ```dataset_parameters.yaml```, which specifies dataset metadata and several paths on your filesystem.

### Preparing sampl