
# Overview
This is a comprehensive pipeline for amplicon-based metagenomics integrating in a Snakemake workflow the [best functions of many tools](#References). It enables performant and reproducibile processing of 16S rRNA or ITS Illumina paired-end reads. The whole process from local .fastq or SRA depository files to generation of basic visualization plots, including quality control plots of intermediate steps, is covered. 
    
# Installation
## Dependencies
* [Mamba](https://github.com/conda-forge/miniforge)
* [Apptainer](https://github.com/apptainer/apptainer/blob/main/INSTALL.md)
## with git
```sh
git clone https://github.com/metagenlab/zAMP.git
pip install -e zAMP
```
# Usage
```sh
$ zamp -h
Usage: zamp [OPTIONS] COMMAND [ARGS]...

  Snakemake pipeline designed for convenient, reproducible and scalable
  amplicon-based metagenomics

  For more options, run: zamp command --help

Options:
  -v, --version  Show the version and exit.
  -h, --help     Show this message and exit.

Commands:
  db        Prepare database files for zAMP
  run       Run zAMP
  citation  Print zAMP and tools citations
```

# Example usage
## Run command
```sh
zamp run -i samples.tsv -db greengenes2 --fw-primer CCTACGGGNGGCWGCAG --rv-primer GACTACHVGGGTATCTAATCC 
```



