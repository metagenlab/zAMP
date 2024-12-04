
# Overview
zAMP is a [Snakemake](https://github.com/snakemake/snakemake) cli, written in [Snaketool](https://github.com/beardymcjohnface/Snaketool), for amplicon metagenomics anlysis. It includes [state-of-the-art tools](zamp/zamp.CITATION) for performant and reproducibile analysis of 16S rRNA or ITS Illumina paired-end reads. 

Starting from local fastq or SRA reads zAMP performs reads QC, ASV inference, taxonomic assignments and basic visualization plots. In addition, zAMP enables training taxonomic classifiers on specific primer-amplified regions in popular databases like [greengenes2](https://greengenes2.ucsd.edu/), [SILVA](https://www.arb-silva.de/) and [UNITE](https://unite.ut.ee/) to increase sensitivity.

## Installation

### Dependencies

* [Mamba](https://github.com/conda-forge/miniforge)
* [Apptainer](https://github.com/apptainer/apptainer/blob/main/INSTALL.md)

### from source
```sh
git clone https://github.com/metagenlab/zAMP.git
pip install -e zAMP
```

## Usage

### Prepare database

**Greengenes2**

* Download

```sh
wget http://ftp.microbio.me/greengenes_release/2022.10/2022.10.backbone.full-length.fna.qza
wget http://ftp.microbio.me/greengenes_release/2022.10/2022.10.backbone.tax.qza
```

* Export qza with [qiime2 export](https://docs.qiime2.org/2024.5/tutorials/exporting/)
```sh
docker run -t -i -v $(pwd):/data quay.io/qiime2/tiny:2024.5 \
qiime tools export \
--input-path 2022.10.backbone.full-length.fna.qza \
--output-path greengenes2
docker run -t -i -v $(pwd):/data quay.io/qiime2/tiny:2024.5 \
qiime tools export \
--input-path 2022.10.backbone.tax.qza --output-path greengenes2
```

* Prepare database 
```sh
zamp db --fasta greengenes2/dna-sequences.fasta \
--taxonomy greengenes2/taxonomy.tsv --name greengenes2 \
--fw-primer CCTACGGGNGGCWGCAG --rv-primer GACTACHVGGGTATCTAATCC \
-o greengenes2
```


### Run with prepared database
```sh
zamp run -i zamp/data/sra-samples.tsv \
-db greengenes2 \
--fw-primer CCTACGGGNGGCWGCAG --rv-primer GACTACHVGGGTATCTAATCC 
```


## Help
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
