import re
import yaml
import os

## When using singularity
if "--use-singularity" in sys.argv:    
    ### Bind the directory of the database to the singularity containers.
    workflow.singularity_args += f' -B {config["tax_DB_path"]}:{config["tax_DB_path"]}'
    #### Load a dictionnary of singularity containers that will be called from each rule
    singularity_envs = yaml.safe_load(open(os.path.join(workflow.basedir,  "envs/singularity/sing_envs.yml"), 'r'))


include: config["assembly_finder_Snakefile"]

## Include the pipeline rules
include: "rules/0_preprocessing/logging.rules"
include: "rules/3_tax_assignment/tax_assign.rules"
include: "rules/5_visualization/QIIME2_import.rules"
include: "rules/PICRUSt2/picrust.rules"
include: "rules/In_silico/insilico_validation.rules"


rule insilico_validation:
    input:
        expand("InSilico/3_classified/{classifier}_{tax_DB}/dna-sequences_tax_assignments.txt" , classifier = config["classifier"], tax_DB = config["tax_DB"]),
        "InSilico/2_denoised/count_table.txt",
        expand("InSilico/3_classified/{classifier}_{tax_DB}/InSilico_compare_tax.tsv", classifier = config["classifier"], tax_DB = config["tax_DB"])
