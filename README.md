# pol-ii-speed
This repository contains a pipeline using total RNA-seq to estimate the changes in RNA Polymerase II
elongation speed, and also the changes in the speed of intron splicing.

For details, please read our paper *Estimating changes in RNA Polymerase II
elongation and intron splicing speeds from total
RNA-seq data*. The pre-print is [available on BioRxiv](https://doi.org/10.1101/2025.08.24.672013).

The code is organized as a Nextflow pipeline. To run it, you will need FASTQ
data of total RNA-seq; please read the documentation below for more details.

The model itself is implemented using PyTorch; the code can be found in the [pol_ii_speed_modeling](pol_ii_speed_modeling) folder 
(which can also be installed as a pip package).

## Running the pipeline
Besides the FASTQ data, you will need to prepare two files to run the pipeline:
 * ```samplesheet.csv```, containing sample annotation.
 * ```dataset_params.yaml```, which specifies dataset metadata and several paths on your filesystem.

### Preparing samplesheet

The file ```samplesheet.csv``` needs to contain the following 4 columns: 
 * *sample*, specifying sample names (the names can be arbitrary).
 * *fastq_1*, name of FASTQ file with reads 1 (compressed by gunzip).
* *fastq_2*, same as above for reads 2.
* *strandedness*, specifying the strandedness (read orientation) of the samples. Currently, only values ```forward``` and ```reverse``` are supported.
  (The most common Illumina RNA-seq protocols use ```reverse``` orientations).

Besides these 4 mandatory columns, the samplesheet can contain arbitrary explanatory variables (e.g., genotype, intervention etc.), that can
be used in the dataset parameter file to specify a design formula (see below).

### Preparing dataset parameters

The file ```dataset_params.yaml``` can be created by copying the [dataset_params_template.yaml](dataset_params_template.yaml) file
and filling it according to the comments. We will now discuss several parameter details:

* **Reference genome and transcriptome** (*gtf_file*, *genome_fasta*, *transcriptome_fasta*, and *gtf_source*): In
the pipeline, we are aligning reads both to the genome and transcriptome. Therefore, the reference
files (*gtf_file*, *genome_fasta*, and *transcriptome_fasta*) need to be all compatible (from the same source and version).

  We are using [STAR](https://github.com/alexdobin/STAR) for aligning reads to the genome; please follow the recommendations
from the STAR documentation for details on how to choose the reference genome files 
(see section *2: Generating genome indexes* of the [manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)).

  Also, since Ensembl and Gencode annotations slightly differ, the used source needs to be stated in the *gtf_source* parameter 
(supported values are ```"ensembl"``` and ```"gencode"```).
* **STAR and Salmon indices**: You can optionally pre-build the STAR index and/or Salmon index yourself and provide the path to it. If no path
is provided (the parameter is left ```null```), the index will be built from the provided genome/transcriptome files.
* **Design formula**: the parameter *design_formula* uses R syntax to generate a design matrix from the formula (see [R documentation](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/formula) for details).
The explanatory variables used in the formula need to be the columns of the ```samplesheet.csv```.
* **Factor reference levels**: We are using a likelihood-ratio test to assess the significance of the LFC parameters (p-value comparing to the null hypothesis that the LFC is 0). However, with our current implementation, the categorical variables (factors)
are only compared to the base level. This means that for factors with 3 or more levels (categories), some comparisons are not made.

  By default, the base level is decided by alphabetical order. The parameter *factor_reference_levels* can optionally contain a path to a CSV file,
which needs to contain columns *factor_name* and *reference*. The entries in *factor_name* are supposed to be the column names of the
factors specified in ```samplesheet.csv```. The corresponding entries in *reference* columns specify the base level.

  For example, assume that ```samplesheet.csv``` contains column *genotype*, with possible values assigned to the samples being *wild_type*, *knockout_1* and *knockout_2*.
By default, the value *knockout_1* will be chosen as the base (reference) level, since it's alphabetically first. The design matrix will then contain variables *knockout_2* and *wild_type* (the reference level, *knockout_1*, is not present in the design matrix to avoid perfect multicollinearity).

  To use *wild_type* as the reference level instead, create a CSV file like this:
  ```
  factor_name,reference
  genotype,wild_type
  ```
  Then, specify the path to this CSV as the *factor_reference_levels* parameter.

* **Stage**: Our pipeline consists of 2 workflows: [pre-processing](./workflows/preprocessing.nf) and [modeling](./workflows/modeling.nf).
The pre-processing contains read alignment and related steps, and is supposed to be
run only once per dataset. The modeling, on the other hand, may be run multiple times, e.g, experimenting with different
design formulas.

  The parameter *stage* specifies whether  both or only one workflow should be run. Possible values are 
```"all"``` (to run both workflows), ```"preprocess"``` and ```"model"```.

### Executing the pipeline

Please install [Nextflow](https://www.nextflow.io/) and [Docker](https://www.docker.com/) on your system.
Then, the pipeline can be run by
```commandline
nextflow run main.nf -params-file dataset_params.yaml
```
See the Nextflow documentation for [details](https://www.nextflow.io/docs/latest/executor.html) on how to run the pipeline on your HPC/cloud system.

**Note on Slurm**: On our HPC system using Slurm, we noticed the following bug: when multiple processes fitting the model are executed on the same cluster node,
their CPU usage somehow collide, resulting in orders of magnitude slower process execution. We suspect that this may be related to the
underlying BLAST setting and/or our cluster set up, and we are currently trying to resolve this issue. In case that you experience similar behavior, please
execute at most one *RunModel* process per cluster node.

## Contact
If you have any questions or experience any problems with the code,
please open an issue on this repository, or reach out at
[jkoubele@uni-koeln.de](mailto:jkoubele@uni-koeln.de).