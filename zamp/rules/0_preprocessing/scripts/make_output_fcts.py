
### Define values of grouping keys
def get_grouping_key(column_of_interest):
    file_list = []
    ## Convert to list if provided as a string in config
    if isinstance(column_of_interest,str):
        column_of_interest=[column_of_interest]
    unique_col=list(set(column_of_interest))
    for i in unique_col:
        combined_values = expand("{column_of_interest}/{column_values}", column_of_interest = i, column_values = list(set(all_samples[i])))
        file_list.extend(combined_values)
    return(file_list)

### Define rarefaction levels (values in config + no rarefaction)
def get_rarefaction_key(rarefaction_values):
    file_list = []
    for i in set(rarefaction_values):
        combined_values = expand("rarefaction_{rarefaction_values}", rarefaction_values = i)
        file_list = file_list + combined_values
    file_list = file_list + ["norarefaction"]

    return(file_list)


# Transform the collapse levels from the plotting taxonomic rank
def get_taxa_collapse_level_key(collapse_level):

    file_list = []
    for i in set(collapse_level):
        if i == "OTU":
            value = 'no_collapse'
        elif i == "Species" :
            value = 'collap_7'
        elif i == "Genus" :
            value = 'collap_6'
        elif i == "Family" :
            value = 'collap_5'
        elif i == "Order" :
            value = 'collap_4'
        elif i == "Class" :
            value = 'collap_3'
        elif i == "Phylum" :
            value = 'collap_2'
        elif i == "Kingdom" :
            value = 'collap_1'
        else :
            raise ValueError("Forbidden value in taxa level for barplots")
        file_list.append(value)
    file_list = file_list + ["no_collapse"]
    return(file_list)



## Set of function to generate list of output from config ##############################################################

### Light output, including QC, count table, consensus sequences and taxonomic assignement
def light_output_list():
    output = []
    output = MultiQC
    output.append(light_output)
    if "DADA2" in config["denoiser"]:
        output.append("DADA2/2_denoised/DADA2_denoising_stats.tsv")
    return(output)

### Basic output, diagnostic plots, KRONA plots and rarefaction curves
def basic_plots_list():
    output = []
    output = light_output_list()
    output.append(basic_plots)
    return(output)

### Complete set of phyloseq, in option including transposed count table and metadata (wide to long)
def phyloseq_output_list():
    output = []
    output = basic_plots_list()
    output.append(phyloseq)
    if config["melted_phyloseq"] == True:
        output.append(phyloseq_melted)
    if config["transposed_tables"] == True:
        output.append(transposed_output)
    return(output)


### Qiime2 outputs, offering interactive visualization of the data
def Qiime2_output_list():
    output = []
    output = basic_plots_list()
    if config["Qiime2_basic_output_visualization"] == True:
        output.append(Qiime2_vis_qzv)
    return(output)

### PICRUSt2 outputs, to predict genomic content based on taxonomic profiles
def PICRUSt2_list():
    output = []
    output = basic_plots_list()
    output.append(picrust2)
    return(output)

### Rule all, the default rule including all default output (not PICRUSt2, since long to compute)
def rule_all_list():
    output = []
    output.append(basic_plots_list())
    output.append(Qiime2_vis_qzv)
    output.append(phyloseq_output_list())
    return(output)
