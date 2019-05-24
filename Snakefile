
## Input/output
### Generate the input sample list
include: "rules/output_definition/making_sample_dataset.rules"
## Include list of callable outputs integrated in fuctions that will be called by the following rules
include: "rules/output_definition/making_output_list_files.rules"

### Defaul rule all. Include all but PICRUSt2 and QIIME2 special outputs
rule all:
    input: rule_all_list()

## Rules to call defined sets of output. For each, we first generate a function calling the combined list of output. Then its is integrated in a easily callable rule
### Only QC of the reads
rule QC:
    input: MultiQC

### Light output, including count table, consensus sequences and taxonomic assignement
rule light_output:
    input: light_output_list()

### Basic output, generated plots numbers, KRONA plots and rarefaction curve
rule basic_output:
    input: basic_plots_list()

### Complete set of phyloseq, in option including transposed count table and metadata (wide to long)
rule phyloseq_output:
    input: phyloseq_output_list()

### Complete set of plots
rule plots_output:
    input: plots_output_list()

### Qiime2 outputs
rule Qiime2_output:
    input: Qiime2_output_list()

### PICRUSt2 outputs
rule PICRUSt2_output:
    input: PICRUSt2_list()


## Include the pipeline rules
include: "rules/common_preprocessing/get_reads.rules"
include: "rules/common_preprocessing/get_sras.rules"
include: "rules/common_preprocessing/QC_raw_reads.rules"
include: "rules/DADA2_ASV/cutadapt_trim.rules"
include: "rules/DADA2_ASV/DADA2_denoising.rules"
include: "rules/DADA2_ASV/QC_DADA2.rules"
include: "rules/vsearch_OTU/PANDAseq_trim_filter_pair.rules"
include: "rules/vsearch_OTU/vsearch_derep_and_clustering.rules"
include: "rules/vsearch_OTU/vsearch_count_tables_and_reformat.rules"
include: "rules/common_tax_assignment/RDP_in_QIIME.rules"
include: "rules/common_tax_tree/tree.rules"
include: "rules/common_visualization/Import_to_QIIME2.rules"
include: "rules/common_visualization/Phyloseq.rules"
include: "rules/common_visualization/General_plotting.rules"
include: "rules/PICRUSt/picrust.rules"
include: "rules/common_visualization/Qiime2_stat_analysis.rules"
include: "rules/common_visualization/Qiime2_plugins.rules"
