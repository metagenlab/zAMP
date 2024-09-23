import re
import yaml
import os
from attrmap import AttrMap
import attrmap.utils as au
from Bio.Seq import Seq
import glob


# Concatenate Snakemake's own log file with the master log file
def copy_log_file():
    files = glob.glob(os.path.join(".snakemake", "log", "*.snakemake.log"))
    if not files:
        return None
    current_log = max(files, key=os.path.getmtime)
    shell("cat " + current_log + " >> " + LOG)


onsuccess:
    copy_log_file()


onerror:
    copy_log_file()


## config file
configfile: os.path.join(workflow.basedir, "config", "config.yaml")


config = AttrMap(config)
config = au.convert_state(config, read_only=True)


include: os.path.join("rules", "0_preprocessing", "directories.smk")
include: os.path.join("rules", "0_preprocessing", "functions.smk")


## Output and log
OUTPUT = config.args.output
LOG = os.path.join(OUTPUT, "zamp.log")

## Input
INSILICO_TAX = config.args.input_tax

## Database args
DBPATH = os.path.dirname(os.path.abspath(config.args.database))
DBNAME = os.listdir(DBPATH)

## Classifier
CLASSIFIER = config.args.classifier

## Assembly finder args
TAXONKIT = config.args.taxonkit
API_KEY = config.args.api_key
LIMIT = config.args.limit
COMPRESSED = config.args.compressed
SOURCE = config.args.source
INCLUDE = config.args.include
TAXON = config.args.taxon
REFERENCE = config.args.reference
ASM_LVL = config.args.assembly_level
ANNOTATED = config.args.annotated
ATYPICAL = config.args.atypical
MAG = config.args.mag
RANK = config.args.rank
NRANK = config.args.nrank

## In-silico PCR tools
PCR_TOOL = config.args.pcr_tool
### In-silico PCR

### Simulate_PCR args
MISMATCH = config.args.mismatch
THREEPRIME = config.args.threeprime

## Cutadapt args
ERRORS = config.args.errors
FW_PRIMER = config.args.fw_primer
RV_PRIMER = config.args.rv_primer
COV = config.args.ampcov
MAXLEN = config.args.maxlen
MINLEN = config.args.minlen

AMPLICON = config.args.amplicon
ADAPTER = ""
if FW_PRIMER and RV_PRIMER:
    FW_PRIMER_COMPL = Seq.reverse_complement(Seq(FW_PRIMER))
    FW_LEN = len(FW_PRIMER)
    RV_LEN = len(RV_PRIMER)
    RV_PRIMER_COMPL = Seq.reverse_complement(Seq(RV_PRIMER))
    FW_COV = round(FW_LEN * COV)
    RV_COV = round(RV_LEN * COV)
    if AMPLICON == "16S":
        ADAPTER = (
            f"{FW_PRIMER};min_overlap={FW_COV}...{RV_PRIMER_COMPL};min_overlap={RV_COV}"
        )
    else:
        ADAPTER = f"{FW_PRIMER};min_overlap={FW_COV};optional...{RV_PRIMER_COMPL};min_overlap={RV_COV}"

## Replace empty tax
REPL_EMPTY = config.args.replace_empty

## Tax assign
CLASSIFIER = config.args.classifier
DENOISER = config.args.denoiser

## When using singularity
if "--use-singularity" in sys.argv:
    ### Mount Database path to singularity containers.
    workflow.deployment_settings.apptainer_args += f" -B {os.path.abspath(config.args.database)}:{os.path.abspath(config.args.database)}"
    #### Load a dictionnary of singularity containers that will be called from each rule
singularity_envs = yaml.safe_load(
    open(os.path.join(workflow.basedir, dir.envs, "singularity", "sing_envs.yml"), "r")
)


## Include the pipeline rules
include: os.path.join("rules", "In_silico", "assembly_finder.rules")
include: os.path.join("rules", "3_tax_assignment", "tax_assign.rules")
include: os.path.join("rules", "5_visualization", "QIIME2_import.rules")
include: os.path.join("rules", "PICRUSt2", "picrust.rules")
include: os.path.join("rules", "In_silico", "insilico_validation.rules")


rule insilico_validation:
    input:
        expand(
            os.path.join(
                dir.out.base,
                "InSilico",
                "3_classified",
                "{classifier}_{tax_DB}",
                "{files}",
            ),
            classifier=CLASSIFIER,
            tax_DB=DBNAME,
            files=["dna-sequences_tax_assignments.txt", "InSilico_compare_tax.tsv"],
        ),
        os.path.join(dir.out.base, "InSilico", "2_denoised", "count_table.tsv"),
