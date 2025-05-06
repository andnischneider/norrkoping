# norrkoping - Analysis code for Schneider et al.

This repository contains R scripts for parts of the preprocessing and analysis done on DNA sequencing and PLFA data for Schneider et al (2024).


## PLFA Analysis

The data needed for this analysis is also included in this repository.

## Sequencing Data analysis

The raw data is stored in the european nucleotide archive (ENA) with accession PRJEB74806, and was preprocessed using dada2 and a modified version of the snakemake workflow available [here](https://github.com/andnischneider/its_workflow). Resulting dada2 output files were further pre-processed using the script 1_Preprocessing.R


For reproductive purposes, the preprocessed data has been made available [here](https://zenodo.org/records/15348074), and can be further processed and analyzed using the R scripts in the results folder of this repository.
