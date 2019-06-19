## Use Docker image with singularity, if "--use-singularity" is used in snakemake command
singularity: "docker://valscherz/amplicon_snakemake_pipeline:0.0.9-1-alpha"

## Input/output
### Include logging rules
include: "rules/0_preprocessing/logging.rules"
## Include rules and functions to set input/output
include: "rules/0_preprocessing/make_output_list_files.rules"


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

### Complete set of plots
rule plots_output:
    input: plots_output_list()
    priority: 45

### PICRUSt2 outputs
rule PICRUSt2_output:
    input: PICRUSt2_list()
    priority: 0



## Include the pipeline rules
include: "rules/0_preprocessing/get_reads.rules"
include: "rules/0_preprocessing/get_sras.rules"
include: "rules/0_preprocessing/QC_raw_reads.rules"
include: "rules/1_2_DADA2_ASVs/1_cutadapt_trim.rules"
include: "rules/1_2_DADA2_ASVs/2_DADA2_denoising.rules"
include: "rules/1_2_vsearch_OTUs/1_PANDAseq_trim_filter_merge.rules"
include: "rules/1_2_vsearch_OTUs/2_vsearch_denoising.rules"
include: "rules/3_tax_assignment/RDP_in_QIIME.rules"
include: "rules/4_post_processing/physeq_processing.rules"
include: "rules/5_visualization/General_plotting.rules"
include: "rules/5_visualization/QIIME2_import.rules"
include: "rules/PICRUSt2/picrust.rules"
