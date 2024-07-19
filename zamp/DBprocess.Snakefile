from Bio.Seq import Seq
import re
import yaml
import os
import attrmap as ap
import attrmap.utils as au
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


# config file
configfile: os.path.join(workflow.basedir, "../", "config", "config.yaml")


config = ap.AttrMap(config)


OUTPUT = config.args.output
LOG = os.path.join(OUTPUT, "zamp.log")

# Read config args
DB_NAME = config.args.name
DBPATH = os.path.join(OUTPUT, DB_NAME)
FASTA = config.args.fasta
TAXONOMY = config.args.taxonomy
CLASSIFIERS = config.args.classifier
RDP_MEM = config.args.rdp_mem
PROCESS = config.args.processing
TAX_COLLAPSE = dict(config.args.tax_collapse)


## Cutadapt args
ERRORS = config.args.max_mismatch
FW_PRIMER = config.args.fw_primer
RV_PRIMER = config.args.rv_primer
COV = config.args.ampcov
MAXLEN = config.args.maxlen
MINLEN = config.args.minlen

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
else:
    ADAPTER = ""


## When using singularity
if "--use-singularity" in sys.argv:
    ### Bind the directory of the database to the singularity containers.
    # workflow.deployment_settings.apptainer_args += f" -B {OUTPUT}:{OUTPUT}"
    ### Bind the workflow directory to the singularity containers.
    workflow.deployment_settings.apptainer_args += (
        f" -B {workflow.basedir}:{workflow.basedir}"
    )
#### Load a dictionnary of singularity containers that will be called from each rule
singularity_envs = yaml.safe_load(
    open(os.path.join(workflow.basedir, "envs/singularity/sing_envs.yml"), "r")
)


## Include rules:
include: "rules/DB_processing/trace_n_log_DB.rules"
include: "rules/DB_processing/format_n_train_classifiers.rules"
include: "rules/DB_processing/RDP_validation.rules"


## Taxonomy database can be skipped by config parameters. Include the right rules based on this parameter.


if PROCESS:

    include: "rules/DB_processing/DB_preprocessing.rules"

else:

    include: "rules/DB_processing/DB_skip_preprocessing.rules"


## Call default DBPATH
rule all:
    input:
        os.path.join(OUTPUT, "DB.hash"),


## Optional DBPATH for RDP training diagnostics
rule RDP_validation:
    input:
        os.path.join(OUTPUT, "RDP/RDP_leave_seq_out_accuracy.txt"),
        os.path.join(OUTPUT, "RDP/RDP_leave_tax_out_accuracy.txt"),
        os.path.join(OUTPUT, "RDP/RDP_cross_validate.txt"),
