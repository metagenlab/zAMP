
## Include set of rules necessary for generating the list of input and output for the pipeline
### Generate the input sample list
include: "rules/common_preprocessing/making_sample_dataset.rules"
### Functions to deal with the output list
include: "rules/output_definition/making_output_list_fcts.rules"
### List of output
include: "rules/output_definition/making_output_list_files.rules"


## Rules to call defined sets of output

### Only QC of the reads
rule QC:
    input: MultiQC

### Light output, including count table, consensus sequences and taxonomic assignement
def light_output():
    output = MultiQC
    output.append(light_output)
    if "DADA2" in config["denoiser"]:
        output.append("DADA2/2_denoised/DADA2_denoising_stats.tsv")
    return(output)

rule light_output:
    input: light_output()


### Basic output, generated plots numbers, KRONA plots and rarefaction curve
def basic_plots():
    output = light_output()
    output.append(basic_plots)
    return(output)

rule basic_output:
    input: basic_plots()


### Complete set of phyloseq, in option including transosed count table and metadata (wide to long)
def phyloseq_output():
    output = basic_plots()
    output.append(phyloseq)
    if config["transposed_tables"] == True:
        output.append(transposed_output)
    return(output)

rule phyloseq_output:
    input: phyloseq_output()

### Complete set of plots
def plots_output():
    output = basic_plots()
    if config["Barplots"] == True:
        output.append(barplots)
    if config["Heatmaps"] == True:
        output.append(barplots)
    if config["Alpha_divs"] == True:
        output.append(barplots)
    if config["Distance_ordinations"] == True:
        output.append(distance_ordinations)
    if config["Constrained_ordinations"] == True:
        output.append(constrained_ordinations)
    if config["Unconstrained_ordinations"] == True:
        output.append(unconstrained_ordinations)
    return(output)

rule plots_output:
    input: plots_output()

### Qiime2 outputs
def Qiime2_outupt():
    output = basic_plots()
    if config["Qiime2_basic_output_visualization"] == True:
        output.append(Qiime2_vis_qzv)
    if config["Volatility"] == True:
        output.append(Qiime2_volatility)
    if config["ANCOM"] == True:
        output.append(Qiime2_ANCOM)
    if "Gradient" in config["Gneiss"]
        output.append(Qiime2_Gneiss_correlation)
    if "Phylogeny" in config["Gneiss"]
        output.append(Qiime2_Gneiss_Phylogeny)
    if "Clustering" in config["Gneiss"]
        output.append(Qiime2_Gneiss_gradient)
    return(output)

rule plots_output:
    input: Qiime2_outupt()

### PICRUSt2 outputs
def PICRUSt2():
    output = basic_plots()
    output.append(picrust2)
    return(output)

rule PICRUSt2_output:
    input: Qiime2_outupt()

### Rule all
def rule_all_list():
    output = basic_plots()
    append.output(phyloseq_output())
    append.output(plots_output())
    append.output(Qiime2_outupt())
    append.output(PICRUSt2())
    return(output)

rule all:
    input:
        rule_all_list(),


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
