from humanfriendly import parse_size
import os
import re
import yaml
from snakemake.utils import min_version
from attrmap import AttrMap
import attrmap.utils as au
import pandas as pd
import glob
from metasnek import fastq_finder
import warnings

min_version("8.10.6")

warnings.filterwarnings(
    "ignore",
    category=UserWarning,
)


# Common functions
def copy_log_file():
    files = glob.glob(os.path.join(".snakemake", "log", "*.snakemake.log"))
    if not files:
        return None
    current_log = max(files, key=os.path.getmtime)
    shell("cat " + current_log + " >> " + LOG)


# DB process functions command functions
def get_mem_mb(mem):
    return int(parse_size(mem) / 10**6)


# Run command functions
def get_raw_reads(wildcards):
    if SAMPLES.loc[wildcards.sample, "paired"]:
        return list(SAMPLES.loc[wildcards.sample, "R1":"R2"])
    else:
        return list(SAMPLES.loc[wildcards.sample, "R1"])


# Insilico command functions
def get_fasta(wildcards):
    checkout = checkpoints.download_assemblies.get(**wildcards).output[0]
    df = pd.read_csv(checkout, sep="\t")
    df["sample"] = [os.path.basename(path).split("_genomic.fna")[0] for path in df.path]
    df.set_index("sample", inplace=True)
    return df.loc[wildcards.sample].path


def list_amplicons(wildcards):
    checkout = checkpoints.download_assemblies.get(**wildcards).output[0]
    df = pd.read_csv(checkout, sep="\t")
    genomes = [os.path.basename(path).split("_genomic.fna")[0] for path in df.path]
    return expand(
        os.path.join(
            dir.out.base, "InSilico", "1a_trimmed_primers", "{sample}_trimmed.fasta"
        ),
        sample=genomes,
    )


def list_samples_counts(wildcards):
    checkout = checkpoints.download_assemblies.get(**wildcards).output[0]
    df = pd.read_csv(checkout, sep="\t")
    genomes = [os.path.basename(path).split("_genomic.fna")[0] for path in df.path]
    return expand(
        os.path.join(
            dir.out.base,
            "InSilico",
            "2_denoised",
            "countSeqs",
            "{sample}_count_table.tsv",
        ),
        sample=genomes,
    )
