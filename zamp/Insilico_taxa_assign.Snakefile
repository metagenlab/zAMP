import re
import yaml
import os
from attrmap import AttrMap
import attrmap.utils as au
from Bio.Seq import Seq
import glob


include: os.path.join("rules", "0_preprocessing", "functions.smk")


onsuccess:
    copy_log_file()


onerror:
    copy_log_file()


## config file
configfile: os.path.join(workflow.basedir, "config", "config.yaml")


config = AttrMap(config)
config = au.convert_state(config, read_only=True)


include: os.path.join("rules", "0_preprocessing", "directories.smk")


## Output and log
OUTPUT = config.args.output
LOG = os.path.join(OUTPUT, "zamp.log")

## Input
INPUT_TAX = config.args.input_tax

## Database args
DBPATH = os.path.abspath(config.args.database)
if config.args.name:
    DBNAME = config.args.name.split(",")

else:
    DBNAME = os.path.basename(os.path.normpath(DBPATH))
    DBPATH = os.path.dirname(DBPATH)

## Classifier
CLASSIFIER = config.args.classifier

## Assembly finder args
AF_ARGS = config.args.af_args
TAXON = config.args.taxon
LIMIT = config.args.limit
ASM_LEVEL = config.args.assembly
ONLY_REF = config.args.only_ref
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

ADAPTER = ""
if FW_PRIMER and RV_PRIMER:
    FW_PRIMER_COMPL = Seq.reverse_complement(Seq(FW_PRIMER))
    FW_LEN = len(FW_PRIMER)
    RV_LEN = len(RV_PRIMER)
    RV_PRIMER_COMPL = Seq.reverse_complement(Seq(RV_PRIMER))
    FW_COV = round(FW_LEN * COV)
    RV_COV = round(RV_LEN * COV)
    ADAPTER = (
        f"{FW_PRIMER};min_overlap={FW_COV}...{RV_PRIMER_COMPL};min_overlap={RV_COV}"
    )

## Replace empty tax
REPL_EMPTY = config.args.replace_empty

## Tax assign
CLASSIFIER = config.args.classifier
DENOISER = config.args.denoiser

## When using singularity
if "--use-singularity" in sys.argv:
    ### Mount Database path to singularity containers.
    workflow.deployment_settings.apptainer_args += f" -B {os.path.abspath(config.args.database)}:{os.path.abspath(config.args.database)}"
    workflow.deployment_settings.apptainer_args += (
        f" -B {workflow.basedir}:{workflow.basedir}"
    )

    #### Load a dictionnary of singularity containers that will be called from each rule
singularity_envs = yaml.safe_load(
    open(os.path.join(workflow.basedir, dir.envs, "singularity", "sing_envs.yml"), "r")
)


## Include the pipeline rules
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
