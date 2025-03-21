from humanfriendly import parse_size
import os
import re
import yaml
from snakemake.utils import min_version
from attrmap import AttrMap
import attrmap.utils as au
import pandas as pd
import numpy as np
import glob
from metasnek import fastq_finder
import warnings
import itertools

min_version("8.10.6")

warnings.filterwarnings(
    "ignore",
    category=UserWarning,
)


wildcard_constraints:
    classifier="[^/._]+",


os.environ["LC_ALL"] = "C.UTF-8"
os.environ["LANG"] = "C.UTF-8"


# Common functions
def copy_log_file():
    files = glob.glob(os.path.join(".snakemake", "log", "*.snakemake.log"))
    if not files:
        return None
    current_log = max(files, key=os.path.getmtime)
    shell("cat " + current_log + " >> " + LOG)


def list_match_dir(path, pattern):
    search_pattern = os.path.join(path, f"*{pattern}*")
    return [
        os.path.basename(d)
        for d in glob.glob(search_pattern)
        if os.path.isdir(d) and not d.startswith(".")
    ]


# DB process functions
def get_mem_mb(mem):
    return int(parse_size(mem) / 10**6)


def get_sintax_tax(row, ranks):
    sintax = []
    for n, tax in enumerate(row["tax"].split(";")):
        sintax.append(f"{ranks[n][0]}:{tax}")
    return f"{row['seq_id']};tax={','.join(sintax)};"


# Run command functions
def get_reads(wildcards):
    if LOCAL:
        if SAMPLES.loc[wildcards.sample, "paired"]:
            return list(SAMPLES.loc[wildcards.sample, "R1":"R2"])
        else:
            return list(SAMPLES.loc[wildcards.sample, "R1"])
    else:
        checkpoint_out = checkpoints.fasterq_dump.get(**wildcards).output[0]
        if SAMPLES.loc[wildcards.sample, "paired"]:
            return [
                os.path.join(dir.out.base, "sra", "{sample}", "{sample}_1.fastq.gz"),
                os.path.join(dir.out.base, "sra", "{sample}", "{sample}_2.fastq.gz"),
            ]
        else:
            return (
                os.path.join(dir.out.base, "sra", "{sample}", "{sample}_1.fastq.gz"),
            )


def sample_list_run_QC(run):
    file_list = []
    samples = list(SAMPLES[SAMPLES["run"] == run].index)
    for sample in samples:
        if SAMPLES.loc[sample, "paired"]:
            suffix = ["R1", "R2"]
        else:
            suffix = ["single"]
        combined_values = expand(
            os.path.join(
                "{qc}", "FastQC", "raw_reads", "{run}", "{sample}_{suffix}_fastqc.zip"
            ),
            qc=dir.out.qc,
            run=run,
            sample=sample,
            suffix=suffix,
        )
        file_list = file_list + combined_values
    return file_list


def sample_list_overall_QC():
    file_list = []

    for run in set(SAMPLES["run"]):
        samples = list(SAMPLES[SAMPLES["run"] == run].index)
        for sample in samples:
            if SAMPLES.loc[sample, "paired"]:
                suffix = ["R1", "R2"]
            else:
                suffix = ["single"]
            combined_values = expand(
                os.path.join(
                    "{qc}",
                    "FastQC",
                    "raw_reads",
                    "{run}",
                    "{sample}_{suffix}_fastqc.zip",
                ),
                qc=dir.out.qc,
                run=run,
                sample=sample,
                suffix=suffix,
            )
            file_list = file_list + combined_values
    return file_list


# Insilico command functions
def get_samples(wildcards):
    if LOCAL:
        SAMPLES = [
            os.path.basename(path).split(f".{suffix}")[0]
            for path in glob.glob(
                os.path.join(os.path.abspath(config.args.input), f"*.{suffix}")
            )
        ]
    else:
        checkout = checkpoints.download_assemblies.get(**wildcards).output[0]
        df = pd.read_csv(checkout, sep="\t")
        SAMPLES = [os.path.basename(path).split(".fna")[0] for path in df.path]
    return SAMPLES


def get_tax_table(wildcards):
    if LOCAL:
        return os.path.abspath(TAX)
    else:
        return [
            os.path.join(dir.out.base, "assembly_finder", "assembly_summary.tsv"),
            os.path.join(dir.out.base, "assembly_finder", "taxonomy.tsv"),
        ]


def get_fasta(wildcards):
    if LOCAL:
        return os.path.join(
            os.path.abspath(config.args.input), f"{wildcards.sample}.{suffix}"
        )
    else:
        checkout = checkpoints.download_assemblies.get(**wildcards).output[0]
        df = pd.read_csv(checkout, sep="\t")
        df["sample"] = [os.path.basename(path).split(".fna")[0] for path in df.path]
        df.set_index("sample", inplace=True)
        return df.loc[wildcards.sample].path


def list_amplicons(wildcards):
    return expand(
        os.path.join(dir.out.base, "in_silico_pcr", "{sample}.fasta"),
        sample=get_samples,
    )


def read_vsearch_outfile(file):
    sample = os.path.basename(file).split("_matches.txt")[0]
    with open(file, "r") as f:
        amps = [line.strip() for line in f.readlines()]
    if not amps:
        amps = ["no_amp"]
    fastas = [sample] * len(amps)
    return pd.DataFrame({"fasta": fastas, "asv_id": amps})


def list_vsearch_matches(wildcards):
    return expand(
        os.path.join(
            dir.out.base,
            "vsearch",
            "matches",
            "{sample}_matches.txt",
        ),
        sample=get_samples,
    )
