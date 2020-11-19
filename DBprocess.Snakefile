
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

## Set logging into the same directory that the DB output by adding logging folder which will be taken in account in the logginf rules
config["logging_folder"] = config["tax_DB_path"] + config["tax_DB_name"] + "/logs/"


## Include rules:
include: "rules/0_preprocessing/scripts/logging.py"
include: "rules/DB_processing/trace_n_log_DB.rules"
include: "rules/DB_processing/format_n_train_classifiers.rules"

## Taxonomy database can be skipped by config parameters
if config["extract_and_merge"] is True:
    include: "rules/DB_processing/DB_preprocessing.rules"
elif config["extract_and_merge"] is False:
    include: "rules/DB_processing/DB_skip_preprocessing.rules"
else:
    raise IOError("'extract_and_merge' must be 'True' or 'False' in config")


include: "rules/DB_processing/RDP_validation.rules"



## Call default output
rule all:
    input:
        config["tax_DB_path"] + config["tax_DB_name"] + "/DB.hash"


## Optionnal output for RDP training diagnostics
rule RDP_validation:
    input:
         config["tax_DB_path"] + config["tax_DB_name"] + "/RDP/RDP_leave_seq_out_accuracy.txt",
         config["tax_DB_path"] + config["tax_DB_name"] + "/RDP/RDP_leave_tax_out_accuracy.txt",
         config["tax_DB_path"] + config["tax_DB_name"] + "/RDP/RDP_cross_validate.txt"