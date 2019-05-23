
## Include set of rules necessary for generating the list of input and output for the pipeline
### Generate the input sample list
include: "rules/common_preprocessing/making_sample_dataset.rules"
### Functions to deal with the output list
include: "rules/common_preprocessing/making_output_list_fcts.rules"
### List of output
include: "rules/common_preprocessing/making_output_list_files.rules"


## Rules to call defined sets of output
rule QC:
    input: MultiQC

rule light_output:
    input: MultiQC + minimal_output

rule minimal_output:
    input: MultiQC + minimal_output + minimal_diagnostic_plots

rule phyloseq_output:
    input: phyloseq_output_list()

rule plots_output:
    input: plots_list()



## Include the needed rules

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
