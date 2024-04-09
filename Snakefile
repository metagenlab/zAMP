
import re
import yaml
import os
from snakemake.utils import min_version

min_version("5.14")



## When using singularity
if "--use-singularity" in sys.argv:    
    ### Bind the directory of the database to the singularity containers.
    workflow.deployment_settings.apptainer_args += f' -B {config["tax_DB_path"]}:{config["tax_DB_path"]}'
    #### Load a dictionnary of singularity containers that will be called from each rule
singularity_envs = yaml.safe_load(open(os.path.join(workflow.basedir,  "envs/singularity/sing_envs.yml"), 'r'))



## Include the pipeline rules
### Sets of script to handle logs, input files and list of output:
include: "rules/0_preprocessing/scripts/logging.py"
include: "rules/0_preprocessing/scripts/make_input_lists.py"
include: "rules/0_preprocessing/scripts/make_output_lists.py"

### Snakemake rules to do the job: 
include: "rules/0_preprocessing/get_inputs.rules"
include: "rules/0_preprocessing/QC_raw_reads.rules"
include: "rules/1_2_DADA2_ASVs/1_cutadapt_trim.rules"
include: "rules/1_2_DADA2_ASVs/2_DADA2_denoising.rules"
include: "rules/1_2_vsearch_OTUs/1_PANDAseq_trim_filter_merge.rules"
include: "rules/1_2_vsearch_OTUs/2_vsearch_denoising.rules"
include: "rules/3_tax_assignment/tax_assign.rules"
include: "rules/4_post_processing/physeq_processing.rules"
include: "rules/4_post_processing/TreeShrink.rules"
include: "rules/5_visualization/General_plotting.rules"
include: "rules/5_visualization/QIIME2_import.rules"
include: "rules/PICRUSt2/picrust.rules"

report: "workflow.rst"

## Rules to call defined sets of output. For each, we first generate a function calling the combined list of output. Then its is integrated in a easily callable rule
### Defaul rule all. Include all but PICRUSt2
rule all:
    input: rule_all_list()
    
### Only QC of the reads
rule QC:
    input: MultiQC
    priority: 50

### Light output, including count table, consensus sequences and taxonomic assignement
rule light_output:
    input: light_output_list()
    priority: 49

### Basic output, generated plots numbers, KRONA plots and rarefaction curve
rule basic_output:
    input: basic_plots_list()
    priority: 48

### Qiime2 outputs
rule Qiime2_output:
    input: Qiime2_output_list()
    priority: 47

### Complete set of phyloseq, in option including transposed count table and metadata (wide to long)
rule phyloseq_output:
    input: phyloseq_output_list()
    priority: 46

### PICRUSt2 outputs
rule PICRUSt2_output:
    input: PICRUSt2_list()
    priority: 0


