import re
import yaml
import os

## When using singularity
if "--use-singularity" in sys.argv:
    ### Bind the directory of the database to the singularity containers.
    workflow.deployment_settings.apptainer_args += f" -B {os.path.abspath(config.args.database)}:{os.path.abspath(config.args.database)}"
    #### Load a dictionnary of singularity containers that will be called from each rule

singularity_envs = yaml.safe_load(
    open(os.path.join(workflow.basedir, "envs/singularity/sing_envs.yml"), "r")
)


## Include Assembly Finder
### Provide default values for a variable required by assembly_finder
config["community_name"] = "genomes"
if config.get("update_ete3", "Not Found") == "Not Found":
    config["update_ete3"] = False


include: "assembly_finder/Snakefile"
## Include the pipeline rules
include: "rules/0_preprocessing/scripts/logging.py"
include: "rules/3_tax_assignment/tax_assign.rules"
include: "rules/5_visualization/QIIME2_import.rules"
include: "rules/PICRUSt2/picrust.rules"
include: "rules/In_silico/insilico_validation.rules"


rule insilico_validation:
    input:
        expand(
            "InSilico/3_classified/{classifier}_{tax_DB}/dna-sequences_tax_assignments.txt",
            classifier=CLASSIFIER,
            tax_DB=DBNAME,
        ),
        "InSilico/2_denoised/count_table.tsv",
        expand(
            "InSilico/3_classified/{classifier}_{tax_DB}/InSilico_compare_tax.tsv",
            classifier=CLASSIFIER,
            tax_DB=DBNAME,
        ),
