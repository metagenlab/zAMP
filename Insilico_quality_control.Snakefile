import re
import yaml
import os


#### Load a dictionnary of singularity containers that will be called from each rule

singularity_envs = yaml.safe_load(
    open(os.path.join(workflow.basedir, "envs/singularity/sing_envs.yml"), "r")
)


## Include the pipeline rules
include: "rules/0_preprocessing/scripts/logging.py"
include: "rules/0_preprocessing/scripts/make_input_lists.py"
include: "rules/0_preprocessing/scripts/make_output_lists.py"
include: "rules/0_preprocessing/get_inputs.rules"
include: "rules/1_2_vsearch_OTUs/1_PANDAseq_trim_filter_merge.rules"
include: "rules/1_2_vsearch_OTUs/2_vsearch_denoising.rules"
include: "rules/5_visualization/QIIME2_import.rules"
include: "rules/In_silico/amplicons_reference.rules"


rule check_amplicon_quality:
    input:
        expand(
            "QualityControl/{denoiser}/compare_quality_table.tsv",
            denoiser=config["denoiser"],
        ),
