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

* Remove intronic reads from FASTQ:
```commandline
sh remove_intronic_reads_from_fastq.sh -i /cellfile/datapublic/jkoubele/data_pol_ii/celegans_mutants/FASTQ -o /cellfile/datapublic/jkoubele/data_pol_ii/celegans_mutants/FASTQ_without_intronic_reads -b /cellfile/datapublic/jkoubele/data_pol_ii/celegans_mutants/intronic_reads
```

* Salmon quant on exonic reads:
```commandline
sh salmon_quant.sh -i /cellfile/datapublic/jkoubele/data_pol_ii/celegans_mutants/FASTQ_without_intronic_reads -o /cellfile/datapublic/jkoubele/data_pol_ii/celegans_mutants/salmon_output -s /cellfile/datapublic/jkoubele/reference_genomes/WBcel235/transcriptome/salmon_index
```

* Compute coverage:
```commandline
sh compute_coverage.sh -i /cellfile/datapublic/jkoubele/data_pol_ii/celegans_mutants/intronic_reads/ -o /cellfile/datapublic/jkoubele/data_pol_ii/celegans_mutants/intron_coverage -g /cellfile/datapublic/jkoubele/reference_genomes/WBcel235/ -f Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.fai
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

* Remove intronic reads from FASTQ:
```commandline
sh remove_intronic_reads_from_fastq.sh -i /cellfile/datapublic/jkoubele/data_pol_ii/drosophila_mutants/FASTQ -o /cellfile/datapublic/jkoubele/data_pol_ii/drosophila_mutants/FASTQ_without_intronic_reads -b /cellfile/datapublic/jkoubele/data_pol_ii/drosophila_mutants/intronic_reads
```

* Salmon quant on exonic reads:
```commandline
sh salmon_quant.sh -i /cellfile/datapublic/jkoubele/data_pol_ii/drosophila_mutants/FASTQ_without_intronic_reads -o /cellfile/datapublic/jkoubele/data_pol_ii/drosophila_mutants/salmon_output -s /cellfile/datapublic/jkoubele/reference_genomes/BDGP6.46/transcriptome/salmon_index
```

* Compute coverage:
```commandline
sh compute_coverage.sh -i /cellfile/datapublic/jkoubele/data_pol_ii/drosophila_mutants/intronic_reads/ -o /cellfile/datapublic/jkoubele/data_pol_ii/drosophila_mutants/intron_coverage -g /cellfile/datapublic/jkoubele/reference_genomes/BDGP6.46/ -f Drosophila_melanogaster.BDGP6.46.dna.toplevel.fa.fai
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

* Remove intronic reads from FASTQ:
```commandline
sh remove_intronic_reads_from_fastq.sh -i /cellfile/datapublic/jkoubele/data_pol_ii/aa5_mouse/FASTQ_trimmed -o /cellfile/datapublic/jkoubele/data_pol_ii/aa5_mouse/FASTQ_without_intronic_reads -b /cellfile/datapublic/jkoubele/data_pol_ii/aa5_mouse/intronic_reads
```

* Salmon quant on exonic reads:
```commandline
sh salmon_quant.sh -i /cellfile/datapublic/jkoubele/data_pol_ii/aa5_mouse/FASTQ_without_intronic_reads -o /cellfile/datapublic/jkoubele/data_pol_ii/aa5_mouse/salmon_output -s /cellfile/datapublic/jkoubele/reference_genomes/GRCm39/transcriptome/salmon_index
```

* Compute coverage:
```commandline
sh compute_coverage.sh -i /cellfile/datapublic/jkoubele/data_pol_ii/aa5_mouse/intronic_reads/ -o /cellfile/datapublic/jkoubele/data_pol_ii/aa5_mouse/intron_coverage -g /cellfile/datapublic/jkoubele/reference_genomes/GRCm39/ -f Mus_musculus.GRCm39.dna.primary_assembly.fa.fai
```

## Mouse myocardium

* Align:
```commandline
sh align.sh -i /cellfile/datapublic/jkoubele/data_pol_ii/mouse_myocardium/FASTQ_trimmed -o /cellfile/datapublic/jkoubele/data_pol_ii/mouse_myocardium/STAR_output -s /cellfile/datapublic/jkoubele/reference_genomes/GRCm39/STAR_index
```

* Extract intronic reads:
```commandline
sh extract_intronic_reads.sh -i /cellfile/datapublic/jkoubele/data_pol_ii/mouse_myocardium/STAR_output -o /cellfile/datapublic/jkoubele/data_pol_ii/mouse_myocardium/intronic_reads -g /cellfile/datapublic/jkoubele/reference_genomes/GRCm39
```

* Remove intronic reads from FASTQ:
```commandline
sh remove_intronic_reads_from_fastq.sh -i /cellfile/datapublic/jkoubele/data_pol_ii/mouse_myocardium/FASTQ_trimmed -o /cellfile/datapublic/jkoubele/data_pol_ii/mouse_myocardium/FASTQ_without_intronic_reads -b /cellfile/datapublic/jkoubele/data_pol_ii/mouse_myocardium/intronic_reads
```

* Salmon quant on exonic reads:
```commandline
sh salmon_quant.sh -i /cellfile/datapublic/jkoubele/data_pol_ii/mouse_myocardium/FASTQ_without_intronic_reads -o /cellfile/datapublic/jkoubele/data_pol_ii/mouse_myocardium/salmon_output -s /cellfile/datapublic/jkoubele/reference_genomes/GRCm39/transcriptome/salmon_index
```

* Compute coverage:
```commandline
sh compute_coverage.sh -i /cellfile/datapublic/jkoubele/data_pol_ii/mouse_myocardium/intronic_reads/ -o /cellfile/datapublic/jkoubele/data_pol_ii/mouse_myocardium/intron_coverage -g /cellfile/datapublic/jkoubele/reference_genomes/GRCm39/ -f Mus_musculus.GRCm39.dna.primary_assembly.fa.fai
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

* Remove intronic reads from FASTQ:
```commandline
sh remove_intronic_reads_from_fastq.sh -i /cellfile/datapublic/jkoubele/data_pol_ii/human_astrocytes/FASTQ -o /cellfile/datapublic/jkoubele/data_pol_ii/human_astrocytes/FASTQ_without_intronic_reads -b /cellfile/datapublic/jkoubele/data_pol_ii/human_astrocytes/intronic_reads
```

* Salmon quant on exonic reads:
```commandline
sh salmon_quant.sh -i /cellfile/datapublic/jkoubele/data_pol_ii/human_astrocytes/FASTQ_without_intronic_reads -o /cellfile/datapublic/jkoubele/data_pol_ii/human_astrocytes/salmon_output -s /cellfile/datapublic/jkoubele/reference_genomes/GRCh38.p14/transcriptome/salmon_index
```

* Compute coverage:
```commandline
sh compute_coverage.sh -i /cellfile/datapublic/jkoubele/data_pol_ii/human_astrocytes/intronic_reads/ -o /cellfile/datapublic/jkoubele/data_pol_ii/human_astrocytes/intron_coverage -g /cellfile/datapublic/jkoubele/reference_genomes/GRCh38.p14/ -f Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai
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

* Remove intronic reads from FASTQ:

```commandline
sh remove_intronic_reads_from_fastq.sh -i /cellfile/datapublic/jkoubele/data_pol_ii/senescent_cell_Apapantonis/FASTQ -o /cellfile/datapublic/jkoubele/data_pol_ii/senescent_cell_Apapantonis/FASTQ_without_intronic_reads -b /cellfile/datapublic/jkoubele/data_pol_ii/senescent_cell_Apapantonis/intronic_reads
```

* Salmon quant on exonic reads:
```commandline
sh salmon_quant.sh -i /cellfile/datapublic/jkoubele/data_pol_ii/senescent_cell_Apapantonis/FASTQ_without_intronic_reads -o /cellfile/datapublic/jkoubele/data_pol_ii/senescent_cell_Apapantonis/salmon_output -s /cellfile/datapublic/jkoubele/reference_genomes/GRCh38.p14/transcriptome/salmon_index
```

* Compute coverage:
```commandline
sh compute_coverage.sh -i /cellfile/datapublic/jkoubele/data_pol_ii/senescent_cell_Apapantonis/intronic_reads/ -o /cellfile/datapublic/jkoubele/data_pol_ii/senescent_cell_Apapantonis/intron_coverage -g /cellfile/datapublic/jkoubele/reference_genomes/GRCh38.p14/ -f Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai
```