# Executed Code

## Celegans mutants

* Align:

```commandline
sh align.sh -i /cellfile/datapublic/jkoubele/data_pol_ii/celegans_mutants/FASTQ -o /cellfile/datapublic/jkoubele/data_pol_ii/celegans_mutants/STAR_output -s /cellfile/datapublic/jkoubele/reference_genomes/WBcel235/STAR_index
```

* Extract intronic reads:

```commandline
sh extract_intronic_reads.sh -i /cellfile/datapublic/jkoubele/data_pol_ii/celegans_mutants/STAR_output -o /cellfile/datapublic/jkoubele/data_pol_ii/celegans_mutants/intronic_reads -g /cellfile/datapublic/jkoubele/reference_genomes/WBcel235
```

## Drosophila mutants

* Align:

```commandline
sh align.sh -i /cellfile/datapublic/jkoubele/data_pol_ii/drosophila_mutants/FASTQ -o /cellfile/datapublic/jkoubele/data_pol_ii/drosophila_mutants/STAR_output -s /cellfile/datapublic/jkoubele/reference_genomes/BDGP6.46/STAR_index
```

* Extract intronic reads:

```commandline
sh extract_intronic_reads.sh -i /cellfile/datapublic/jkoubele/data_pol_ii/drosophila_mutants/STAR_output -o /cellfile/datapublic/jkoubele/data_pol_ii/drosophila_mutants/intronic_reads -g /cellfile/datapublic/jkoubele/reference_genomes/BDGP6.46
```

## AA5 mouse

* Align:

```commandline
sh align.sh -i /cellfile/datapublic/jkoubele/data_pol_ii/aa5_mouse/FASTQ_trimmed -o /cellfile/datapublic/jkoubele/data_pol_ii/aa5_mouse/STAR_output -s /cellfile/datapublic/jkoubele/reference_genomes/GRCm39/STAR_index
```

* Extract intronic reads:

```commandline
sh extract_intronic_reads.sh -i /cellfile/datapublic/jkoubele/data_pol_ii/aa5_mouse/STAR_output -o /cellfile/datapublic/jkoubele/data_pol_ii/aa5_mouse/intronic_reads -g /cellfile/datapublic/jkoubele/reference_genomes/GRCm39
```

## Human Astrocytes

* Align:

```commandline
sh align.sh -i /cellfile/datapublic/jkoubele/data_pol_ii/human_astrocytes/FASTQ -o /cellfile/datapublic/jkoubele/data_pol_ii/human_astrocytes/STAR_output -s /cellfile/datapublic/jkoubele/reference_genomes/GRCh38.p14/STAR_index
```

* Extract intronic reads:

```commandline
sh extract_intronic_reads.sh -i /cellfile/datapublic/jkoubele/data_pol_ii/human_astrocytes/STAR_output -o /cellfile/datapublic/jkoubele/data_pol_ii/human_astrocytes/intronic_reads -g /cellfile/datapublic/jkoubele/reference_genomes/GRCh38.p14
```

## Human senescent cells (Apapantonis data AKI46)

* Align:

```commandline
sh align.sh -i /cellfile/datapublic/jkoubele/data_pol_ii/senescent_cell_Apapantonis/FASTQ -o /cellfile/datapublic/jkoubele/data_pol_ii/senescent_cell_Apapantonis/STAR_output -s /cellfile/datapublic/jkoubele/reference_genomes/GRCh38.p14/STAR_index
```

* Extract intronic reads:

```commandline
sh extract_intronic_reads.sh -i /cellfile/datapublic/jkoubele/data_pol_ii/senescent_cell_Apapantonis/STAR_output -o /cellfile/datapublic/jkoubele/data_pol_ii/senescent_cell_Apapantonis/intronic_reads -g /cellfile/datapublic/jkoubele/reference_genomes/GRCh38.p14
```