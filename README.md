# pol-ii-speed
This repository contain a pipeline using total RNA-seq to estimate the changes of RNA Polymerase II
elongation speed, and also the changes in the speed of intron splicing.

For details, please read our paper *Estimating changes in RNA Polymerase II
elongation and intron splicing speeds from total
RNA-seq data*. The pre-print should be on a BioRxiv; the PDF is also available in this repository [here](./paper.pdf).

The code is organized as a Nextflow pipeline. To run it, you will need FASTQ
data of total RNA-seq; please read the documentation below for more details.

The model itself is implemented using PyTorch; the code can be found in the [pol_ii_speed_modeling](pol_ii_speed_modeling) folder 
(which can be also installed as a pip package).

## Running the pipeline
Besides the FASTQ data, you will need to prepare two files to run the pipeline:
 * ```samplesheet.csv```, containing sample annotation.
 * ```dataset_params.yaml```, which specifies dataset metadata and several paths on your filesystem.

### Preparing samplesheet

The file ```samplesheet.csv``` needs to containg following 4 columns: 
 * *sample*, specifying sample names (the names can be arbitrary).
 * *fastq_1*, name of FASTQ file with reads 1 (compressed by gunzip).
* *fastq_2*, same as above for reads 2.
* *strandedness*, specifying the strandedness (read orientation) of the samples. Currently, only values ```forward``` and ```reverse``` are supported.
  (The most common Illumina RNA-seq protocols use ```reverse``` orientations).

Besides these 4 mandatory columns, the samplesheet can contain arbitrary explanatory variables (e.g., genotype, intervention etc.), that can
be used in the dataset parameter file to specify a design formula (see below).

### Preparing dataset parameters

The file ```dataset_params.yaml``` can be created by copying the [dataset_params_template.yaml](dataset_params_template.yaml) file
and filling it according to the comments.

### Executing the pipeline

Please install [Nextflow](https://www.nextflow.io/) and [Docker](https://www.docker.com/) on your system.
Then, the pipeline can be run by
```commandline
nextflow run main.nf -params-file dataset_params.yaml
```