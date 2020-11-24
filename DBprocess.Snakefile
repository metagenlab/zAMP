
import re
import yaml
import os

## When using singularity
if "--use-singularity" in sys.argv:    
    ### Bind the directory of the database to the singularity containers.
    workflow.singularity_args += f' -B {config["tax_DB_path"]}:{config["tax_DB_path"]}'
    ### Bind the workflow directory to the singularity containers.
    workflow.singularity_args += f' -B {workflow.basedir}:{workflow.basedir}'
#### Load a dictionnary of singularity containers that will be called from each rule
singularity_envs = yaml.safe_load(open(os.path.join(workflow.basedir,  "envs/singularity/sing_envs.yml"), 'r'))


## Format output path
path_to_DB = config["tax_DB_path"]
if not path_to_DB.endswith("/"):
    path_to_DB = path_to_DB + "/"

processed_DB_dir = config["tax_DB_name"]
if not processed_DB_dir.endswith("/"):
    processed_DB_dir = processed_DB_dir + "/"

processed_DB_path = path_to_DB + processed_DB_dir


## Set logging into the same directory that the DB output by adding logging folder which will be taken in account in the logginf rules
config["logging_folder"] = processed_DB_path + "/logs/"


## Include rules:
include: "rules/0_preprocessing/scripts/logging.py"
include: "rules/DB_processing/trace_n_log_DB.rules"
include: "rules/DB_processing/format_n_train_classifiers.rules"
include: "rules/DB_processing/RDP_validation.rules"

## Taxonomy database can be skipped by config parameters. Include the right rules based on this parameter.
if config["extract_and_merge"] is True:
    include: "rules/DB_processing/DB_preprocessing.rules"
elif config["extract_and_merge"] is False:
    include: "rules/DB_processing/DB_skip_preprocessing.rules"
else:
    raise IOError("'extract_and_merge' must be 'True' or 'False' in config")


## Call default output
rule all:
    input:
        processed_DB_path + "DB.hash"


## Optional output for RDP training diagnostics
rule RDP_validation:
    input:
         processed_DB_path + "RDP/RDP_leave_seq_out_accuracy.txt",
         processed_DB_path + "RDP/RDP_leave_tax_out_accuracy.txt",
         processed_DB_path + "RDP/RDP_cross_validate.txt"