include: "rules/common_preprocessing/making_sample_dataset.rules"

## Function to modulate keys


def get_alpha_diversity_files(config):

    import itertools

    # list of factors
    denoiser = config["denoiser"]
    tax_DB = config["tax_DB"]
    filter_tax_rank = config["filter_tax_rank"]
    rarefaction_value = ["rarefaction_%s" %i for i in config["rarefaction_value"]] + ["norarefaction"]
    filter_lineage = config["filter_lineage"]
    filter_column_value = config["filter_column_value"]
    filter_meta_column = config["filter_meta_column"]

    list_of_lists = [denoiser, tax_DB, rarefaction_value, filter_tax_rank,filter_lineage, filter_column_value, filter_meta_column]
    complete_coombination_list = list(itertools.product(*list_of_lists))

    # colums from which factors have to be extracted
    grouping_column = list(set(all_samples[config["grouping_column"]]))

    file_list = []

    template = '%s/5_visualization/rdp/%s/%s/alpha_diversity/%s_%s_taxfilt_%s_in_%s/%s/%s_alpha_divesity.png'
    for column in grouping_column:
        grouping_factor_list = list(set(all_samples[column]))
        for combination in complete_coombination_list:
            for grouping_factor in grouping_factor_list:
                denoiser = combination[0]
                tax_DB = combination[1]
                rarefaction_value = combination[2]
                filter_tax_rank = combination[3]
                filter_lineage = combination[4]
                filter_column_value = combination[5]
                filter_meta_column = combination[6]
                #print(combination)

                file_name = template % (denoiser,
                                        tax_DB,
                                        rarefaction_value,
                                        filter_tax_rank,
                                        filter_lineage,
                                        filter_column_value,
                                        filter_meta_column,
                                        column,
                                        grouping_factor)
                #print(file_name)
                file_list.append(file_name)


    #template = '{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/alpha_diversity/{filter_tax_rank}_{filter_lineage}_taxfilt_{filter_column_value}_in_{filter_meta_column}/{grouping_column}/{grouping_column_value}_alpha_divesity.png'
    return file_list


def get_grouping_key(column_of_interest):

    filtered_all_samples = all_samples.loc[all_samples[config['filter_meta_column']]  == config['filter_column_value']]

    file_list = []

    for i in set(column_of_interest):
        combined_values = expand("{column_of_interest}/{column_values}", column_of_interest = i, column_values = list(set(filtered_all_samples[i])))
        file_list = file_list + combined_values
    return(file_list)


def get_rarefaction_key(rarefaction_values):
    file_list = []



    for i in set(rarefaction_values):
        combined_values = expand("rarefaction_{rarefaction_values}", rarefaction_values = i)
        file_list = file_list + combined_values
    file_list = file_list + ["norarefaction"]

    return(file_list)


def get_taxa_collapse_level_key(collapse_level):
    file_list = []

    for i in set(collapse_level):
        combined_values = expand("collap_{collapse_level}", collapse_level = i)
        file_list = file_list + combined_values
    file_list = file_list + ["no_collapse"]

    return(file_list)


def get_filtering_key(filtering):
    file_list = []
    for i in set(filtering):
        print(i)
        if i == "nofiltering" :
            filt = ["nofiltering_0"]
        elif i == "absolute" :
            filt =  expand("{filter}_{filtering_value}" , filter = i, filtering_value = config["absolute_filtering_value"])
        elif i == "relative" :
            filt =  expand("{filter}_{filtering_value}" , filter = i, filtering_value = config["relative_filtering_value"])
        else :
            raise ValueError("Forbidden value for filtering type")

        file_list = file_list + filt

    return(file_list)



## A function do define the ouput based on the content of the config file:
def get_final_output(config):

    ## Base output
    lst=[

    ## QC
    "QC/multiqc_raw_reads_report.html",
    #"QC/multiqc_DADA2_filtered_reads_report.html",

    ## Generic output
    expand("{denoiser}/2_denoised/dna-sequences.fasta", denoiser = config["denoiser"]),
    expand("{denoiser}/2_denoised/count-table.qzv", denoiser = config["denoiser"]),
    expand("{denoiser}/2_denoised/rep-seqs.qzv", denoiser = config["denoiser"]),
    expand("{denoiser}/3_classified/rdp/{tax_DB}/dna-sequences_tax_assignments.qzv", denoiser = config["denoiser"], tax_DB = config["tax_DB"]),
    expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/reads/reads_plot_with_filtered.png",
           denoiser = config["denoiser"],
           tax_DB = config["tax_DB"]),


    ### Rarefaction curve
    expand("{denoiser}/5_visualization/rdp/{tax_DB}/{rarefaction_value}/rarefaction_curve.png",
        denoiser = config["denoiser"],
        tax_DB = config["tax_DB"],
        rarefaction_value = get_rarefaction_key(config["rarefaction_value"])),

    ### KRONA
    expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/KRONA/{grouping_key}.html",
        denoiser = config["denoiser"],
        tax_DB = config["tax_DB"],
        rarefaction_value = get_rarefaction_key(config["rarefaction_value"]),
        grouping_key = get_grouping_key(config["grouping_column"])),

    ### Alpha diversity
    expand("{denoiser}/5_visualization/rdp/{tax_DB}/{rarefaction_value}/alpha_diversity/{filter_tax_rank}_{filter_lineage}_taxfilt_{filter_column_value}_in_{filter_meta_column}/{grouping_key}_alpha_divesity.png",
        denoiser = config["denoiser"],
        tax_DB = config["tax_DB"],
        rarefaction_value = get_rarefaction_key(config["rarefaction_value"]),
        filter_tax_rank = config["filter_tax_rank"],
        filter_lineage = config["filter_lineage"],
        filter_column_value = config["filter_column_value"],
        filter_meta_column = config["filter_meta_column"],
        grouping_key = get_grouping_key(config["grouping_column"])),


    ### Ordination
        #### Distance-based
        expand("{denoiser}/5_visualization/rdp/{tax_DB}/{rarefaction_value}/ordination/{filter_tax_rank}_{filter_lineage}_taxfilt_{filter_column_value}_in_{filter_meta_column}/distance_based/{ordination_method_distance_based}/{grouping_key}_d_{ordination_distance}.png",
            denoiser = config["denoiser"],
            tax_DB = config["tax_DB"],
            rarefaction_value = get_rarefaction_key(config["rarefaction_value"]),
            filter_tax_rank = config["filter_tax_rank"],
            filter_lineage = config["filter_lineage"],
            filter_column_value = config["filter_column_value"],
            filter_meta_column = config["filter_meta_column"],
            ordination_method_distance_based = config["ordination_method_distance_based"],
            grouping_key = get_grouping_key(config["grouping_column"]),
            ordination_distance = config["ordination_distance"]),

        #### Unconstrained
        expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/ordination/{filter_tax_rank}_{filter_lineage}_taxfilt_{filter_column_value}_in_{filter_meta_column}/unconstrained/{ordination_method_unconstrained}/{grouping_key}.png",
            denoiser = config["denoiser"],
            tax_DB = config["tax_DB"],
            filter_tax_rank = config["filter_tax_rank"],
            filter_lineage = config["filter_lineage"],
            filter_column_value = config["filter_column_value"],
            filter_meta_column = config["filter_meta_column"],
            ordination_method_unconstrained = config["ordination_method_unconstrained"],
            grouping_key = get_grouping_key(config["grouping_column"])),

        #### Constrained
        expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/ordination/{filter_tax_rank}_{filter_lineage}_taxfilt_{filter_column_value}_in_{filter_meta_column}/constrained/{ordination_method_constrained}/{grouping_key}_f_{ordination_factor}.png",
            denoiser = config["denoiser"],
            tax_DB = config["tax_DB"],
            filter_tax_rank = config["filter_tax_rank"],
            filter_lineage = config["filter_lineage"],
            filter_column_value = config["filter_column_value"],
            filter_meta_column = config["filter_meta_column"],
            ordination_method_constrained = config["ordination_method_constrained"],
            grouping_key = get_grouping_key(config["grouping_column"]),
            ordination_factor = config["ordination_factor"]),

        ### Melted Phyloseq
        expand("{denoiser}/5_visualization/rdp/{tax_DB}/{rarefaction_value}/physeq/{collapse_key}/2_filter_samples/{filter_tax_rank}_{filter_lineage}_taxfilt_{filter_column_value}_in_{filter_meta_column}_samples_melted.tsv",
            denoiser = config["denoiser"],
            tax_DB = config["tax_DB"],
            collapse_key = get_taxa_collapse_level_key(config["collapse_level"]),
            rarefaction_value = get_rarefaction_key(config["rarefaction_value"]),
            filter_tax_rank = config["filter_tax_rank"],
            filter_lineage = config["filter_lineage"],
            filter_column_value = config["filter_column_value"],
            filter_meta_column = config["filter_meta_column"]),

        ### Melted Phyloseq - percent transformed
        expand("{denoiser}/5_visualization/rdp/{tax_DB}/{rarefaction_value}/physeq/{collapse_key}/2_filter_samples/{filter_tax_rank}_{filter_lineage}_taxfilt_{filter_column_value}_in_{filter_meta_column}_samples_trfs_melted.tsv",
            denoiser = config["denoiser"],
            tax_DB = config["tax_DB"],
            collapse_key = get_taxa_collapse_level_key(config["collapse_level"]),
            rarefaction_value = get_rarefaction_key(config["rarefaction_value"]),
            filter_tax_rank = config["filter_tax_rank"],
            filter_lineage = config["filter_lineage"],
            filter_column_value = config["filter_column_value"],
            filter_meta_column = config["filter_meta_column"]),

        ### Count-table transposed format
        expand("{denoiser}/5_visualization/rdp/{tax_DB}/{rarefaction_value}/physeq/{collapse_key}/3_filter_features/{filter_tax_rank}_{filter_lineage}_taxfilt_{filter_column_value}_in_{filter_meta_column}_features_trfs_export/count_table_transposed.txt",
            denoiser = config["denoiser"],
            tax_DB = config["tax_DB"],
            collapse_key = get_taxa_collapse_level_key(config["collapse_level"]),
            rarefaction_value = get_rarefaction_key(config["rarefaction_value"]),
            filter_tax_rank = config["filter_tax_rank"],
            filter_lineage = config["filter_lineage"],
            filter_column_value = config["filter_column_value"],
            filter_meta_column = config["filter_meta_column"]),

        ### Count-table transposed format
        expand("{denoiser}/5_visualization/rdp/{tax_DB}/{rarefaction_value}/physeq/{collapse_key}/2_filter_samples/{filter_tax_rank}_{filter_lineage}_taxfilt_{filter_column_value}_in_{filter_meta_column}_samples_export/count_table_transposed.txt",
            denoiser = config["denoiser"],
            tax_DB = config["tax_DB"],
            collapse_key = get_taxa_collapse_level_key(config["collapse_level"]),
            rarefaction_value = get_rarefaction_key(config["rarefaction_value"]),
            filter_tax_rank = config["filter_tax_rank"],
            filter_lineage = config["filter_lineage"],
            filter_column_value = config["filter_column_value"],
            filter_meta_column = config["filter_meta_column"]),

        ### Count-table transposed format - filtered
        expand("{denoiser}/5_visualization/rdp/{tax_DB}/{rarefaction_value}/physeq/{collapse_key}/2_filter_samples/{filter_tax_rank}_{filter_lineage}_taxfilt_{filter_column_value}_in_{filter_meta_column}_samples_trfs_export/count_table_transposed.txt",
            denoiser = config["denoiser"],
            tax_DB = config["tax_DB"],
            collapse_key = get_taxa_collapse_level_key(config["collapse_level"]),
            rarefaction_value = get_rarefaction_key(config["rarefaction_value"]),
            filter_tax_rank = config["filter_tax_rank"],
            filter_lineage = config["filter_lineage"],
            filter_column_value = config["filter_column_value"],
            filter_meta_column = config["filter_meta_column"]),

        ### Count-table transposed format - non filtered
        expand("{denoiser}/5_visualization/rdp/{tax_DB}/{rarefaction_value}/physeq/no_collapse/base_export/count_table_transposed.txt",
            denoiser = config["denoiser"],
            tax_DB = config["tax_DB"],
            collapse_key = get_taxa_collapse_level_key(config["collapse_level"]),
            rarefaction_value = get_rarefaction_key(config["rarefaction_value"]),
            filter_tax_rank = config["filter_tax_rank"],
            filter_lineage = config["filter_lineage"],
            filter_column_value = config["filter_column_value"],
            filter_meta_column = config["filter_meta_column"])
    ]

    ## Conditional output

    if config["Barplot"] == True:
        lst.append(
        ### Barplot
        expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/barplot/{filter_tax_rank}_{filter_lineage}_taxfilt_{filter_column_value}_in_{filter_meta_column}/{relative_or_absolute_plot}/{grouping_key}_{filtering_key}_{plotting_tax_ranks}_barplot.png",
        denoiser = config["denoiser"],
        tax_DB = config["tax_DB"],
        rarefaction_value = get_rarefaction_key(config["rarefaction_value"]),
        filter_tax_rank = config["filter_tax_rank"],
        filter_lineage = config["filter_lineage"],
        filter_column_value = config["filter_column_value"],
        filter_meta_column = config["filter_meta_column"],
        relative_or_absolute_plot = config["relative_or_absolute_baxplot"],
        grouping_key = get_grouping_key(config["grouping_column"]),
        filtering_key = get_filtering_key(config["relative_or_absolute_filtering"]),
        plotting_tax_ranks = config["plotting_tax_ranks"]))

    if config["Heatmap"] == True:
        lst.append(
        ### Heatmap
        expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/heatmaps/{filter_tax_rank}_{filter_lineage}_taxfilt_{filter_column_value}_in_{filter_meta_column}/{relative_or_absolute_plot}/{grouping_key}_{filtering_key}_{plotting_tax_ranks}_heatmap.png",
        denoiser = config["denoiser"],
        tax_DB = config["tax_DB"],
        rarefaction_value = get_rarefaction_key(config["rarefaction_value"]),
        filter_tax_rank = config["filter_tax_rank"],
        filter_lineage = config["filter_lineage"],
        filter_column_value = config["filter_column_value"],
        filter_meta_column = config["filter_meta_column"],
        relative_or_absolute_plot = config["relative_or_absolute_baxplot"],
        grouping_key = get_grouping_key(config["grouping_column"]),
        filtering_key = get_filtering_key(config["relative_or_absolute_filtering"]),
        plotting_tax_ranks = config["plotting_tax_ranks"]))


    if config["Volatility"] == True:
       lst.append(
        ## Volatility viz
        expand("{denoiser}/5_visualization/rdp/{tax_DB}/{rarefaction_value}/volatility/{collapse_key}/2_filter_samples/{filter_tax_rank}_{filter_lineage}_taxfilt_{filter_column_value}_in_{filter_meta_column}_samples/volatility.qzv",
            denoiser = config["denoiser"],
            tax_DB = config["tax_DB"],
            collapse_key = get_taxa_collapse_level_key(config["collapse_level"]),
            rarefaction_value = get_rarefaction_key(config["rarefaction_value"]),
            filter_tax_rank = config["filter_tax_rank"],
            filter_lineage = config["filter_lineage"],
            filter_column_value = config["filter_column_value"],
            filter_meta_column = config["filter_meta_column"]))

        ## Volatility viz {tool}/5_visualization/{classifier}/{db_taxonomy}/{raref_or_not}/Qiime2/{collapsed_or_not}/{prefix1}/{prefix2}_export/volatility.qzv
       lst.append(
        expand("{denoiser}/5_visualization/rdp/{tax_DB}/{rarefaction_value}/volatility/{collapse_key}/2_filter_samples/{filter_tax_rank}_{filter_lineage}_taxfilt_{filter_column_value}_in_{filter_meta_column}_samples/feature-volatility_filtered-table.qza",
            denoiser = config["denoiser"],
            tax_DB = config["tax_DB"],
            collapse_key = get_taxa_collapse_level_key(config["collapse_level"]),
            rarefaction_value = get_rarefaction_key(config["rarefaction_value"]),
            filter_tax_rank = config["filter_tax_rank"],
            filter_lineage = config["filter_lineage"],
            filter_column_value = config["filter_column_value"],
            filter_meta_column = config["filter_meta_column"]))

    ## Statistical analyses
    if config["ANCOM"] == True:
        lst.append(
        ## ANCOM
            expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/diff_abundance/{collapse_key}/ANCOM/{filter_tax_rank}_{filter_lineage}_taxfilt_{filter_column_value}_in_{filter_meta_column}_f_{tested_factor}.qzv",
                denoiser = config["denoiser"],
                tax_DB = config["tax_DB"],
                collapse_key = get_taxa_collapse_level_key(config["collapse_level"]),
                filter_tax_rank = config["filter_tax_rank"],
                filter_lineage = config["filter_lineage"],
                filter_column_value = config["filter_column_value"],
                filter_meta_column = config["filter_meta_column"],
                tested_factor= config["ANCOM_factors"]))

    else :
        print("ANCOM not 'True', will not be in output")

    if "Correlation" in config["Gneiss"] :
        lst.append(
        ### Gneiss
        #### Correlation based
        ##### Regression
        expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/diff_abundance/{collapse_key}/Gneiss/correlation/{filter_tax_rank}_{filter_lineage}_taxfilt_{filter_column_value}_in_{filter_meta_column}/hier_correlation_regression.qzv",
            denoiser = config["denoiser"],
            tax_DB = config["tax_DB"],
            collapse_key = get_taxa_collapse_level_key(config["collapse_level"]),
            filter_tax_rank = config["filter_tax_rank"],
            filter_lineage = config["filter_lineage"],
            filter_column_value = config["filter_column_value"],
            filter_meta_column = config["filter_meta_column"])
        )


        ##### Heatmaps
        lst.append(
        expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/diff_abundance/{collapse_key}/Gneiss/correlation/{filter_tax_rank}_{filter_lineage}_taxfilt_{filter_column_value}_in_{filter_meta_column}/hier_correlation_heatmap_{tested_factor}.qzv",
            denoiser = config["denoiser"],
            tax_DB = config["tax_DB"],
            collapse_key = get_taxa_collapse_level_key(config["collapse_level"]),
            filter_tax_rank = config["filter_tax_rank"],
            filter_lineage = config["filter_lineage"],
            filter_column_value = config["filter_column_value"],
            filter_meta_column = config["filter_meta_column"],
            tested_factor= config["ANCOM_factors"])
        )


        ##### Balances
        lst.append(
        expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/diff_abundance/{collapse_key}/Gneiss/correlation/{filter_tax_rank}_{filter_lineage}_taxfilt_{filter_column_value}_in_{filter_meta_column}/hier_correlation_y_{y_balances}_f_{tested_factor}.qzv",
            denoiser = config["denoiser"],
            tax_DB = config["tax_DB"],
            collapse_key = get_taxa_collapse_level_key(config["collapse_level"]),
            filter_tax_rank = config["filter_tax_rank"],
            filter_lineage = config["filter_lineage"],
            filter_column_value = config["filter_column_value"],
            filter_meta_column = config["filter_meta_column"],
            tested_factor= config["ANCOM_factors"],
            y_balances = list(range(1, 9)))
        )



    if "Phylogeny" in config["Gneiss"]:
        ### Phylogeny based
        ##### Regression
        lst.append(
        expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/diff_abundance/{collapse_key}/Gneiss/phylogeny/{filter_tax_rank}_{filter_lineage}_taxfilt_{filter_column_value}_in_{filter_meta_column}/phyl_phylogenetic_regression.qzv",
            denoiser = config["denoiser"],
            tax_DB = config["tax_DB"],
            collapse_key = get_taxa_collapse_level_key(config["collapse_level"]),
            filter_tax_rank = config["filter_tax_rank"],
            filter_lineage = config["filter_lineage"],
            filter_column_value = config["filter_column_value"],
            filter_meta_column = config["filter_meta_column"])
        )


        ##### Heatmaps
        lst.append(
        expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/diff_abundance/{collapse_key}/Gneiss/phylogeny/{filter_tax_rank}_{filter_lineage}_taxfilt_{filter_column_value}_in_{filter_meta_column}/phyl_phylogenetic_heatmap_{tested_factor}.qzv",
            denoiser = config["denoiser"],
            tax_DB = config["tax_DB"],
            collapse_key = get_taxa_collapse_level_key(config["collapse_level"]),
            y_balances = list(range(1, 9)),
            filter_tax_rank = config["filter_tax_rank"],
            filter_lineage = config["filter_lineage"],
            filter_column_value = config["filter_column_value"],
            filter_meta_column = config["filter_meta_column"],
            tested_factor = config["ANCOM_factors"])
        )


        ##### Balances
        lst.append(
        expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/diff_abundance/{collapse_key}/Gneiss/phylogeny/{filter_tax_rank}_{filter_lineage}_taxfilt_{filter_column_value}_in_{filter_meta_column}/phyl_phylogenetic_y_{y_balances}_f_{tested_factor}.qzv",
            denoiser = config["denoiser"],
            tax_DB = config["tax_DB"],
            collapse_key = get_taxa_collapse_level_key(config["collapse_level"]),
            filter_tax_rank = config["filter_tax_rank"],
            filter_lineage = config["filter_lineage"],
            filter_column_value = config["filter_column_value"],
            filter_meta_column = config["filter_meta_column"],
            tested_factor= config["ANCOM_factors"],
            y_balances = list(range(1, 9)))
        )

    if "Gradient" in config["Gneiss"]:
        ### Gradient based
        ##### Regression
        lst.append(
        expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/diff_abundance/{collapse_key}/Gneiss/gradient/{filter_tax_rank}_{filter_lineage}_taxfilt_{filter_column_value}_in_{filter_meta_column}/hier_{tested_factor}_regression.qzv",
            denoiser = config["denoiser"],
            tax_DB = config["tax_DB"],
            collapse_key = get_taxa_collapse_level_key(config["collapse_level"]),
            filter_tax_rank = config["filter_tax_rank"],
            filter_lineage = config["filter_lineage"],
            filter_column_value = config["filter_column_value"],
            filter_meta_column = config["filter_meta_column"],
            tested_factor = config["Gneiss_gradient_clustering"])
        )


        ##### Heatmaps
        lst.append(
        expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/diff_abundance/{collapse_key}/Gneiss/gradient/{filter_tax_rank}_{filter_lineage}_taxfilt_{filter_column_value}_in_{filter_meta_column}/hier_{tested_factor}_heatmap_{tested_factor}.qzv",
            denoiser = config["denoiser"],
            tax_DB = config["tax_DB"],
            collapse_key = get_taxa_collapse_level_key(config["collapse_level"]),
            y_balances = list(range(1, 9)),
            filter_tax_rank = config["filter_tax_rank"],
            filter_lineage = config["filter_lineage"],
            filter_column_value = config["filter_column_value"],
            filter_meta_column = config["filter_meta_column"],
            tested_factor = config["Gneiss_gradient_clustering"])
        )


        ##### Balances
        lst.append(
        expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/diff_abundance/{collapse_key}/Gneiss/gradient/{filter_tax_rank}_{filter_lineage}_taxfilt_{filter_column_value}_in_{filter_meta_column}/hier_{tested_factor}_y_{y_balances}_f_{tested_factor}.qzv",
            denoiser = config["denoiser"],
            tax_DB = config["tax_DB"],
            collapse_key = get_taxa_collapse_level_key(config["collapse_level"]),
            filter_tax_rank = config["filter_tax_rank"],
            filter_lineage = config["filter_lineage"],
            filter_column_value = config["filter_column_value"],
            filter_meta_column = config["filter_meta_column"],
            tested_factor = config["Gneiss_gradient_clustering"],
            y_balances = list(range(1, 9)))
        )

    else :
        print("Gneiss not 'True', will not be in output")



    ###Gneiss gradient based to be implemented


    if config["PICRUST"] == True:
        lst.append(
        ### PIRCUST 2
            expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/picrust/",
                    denoiser = config["denoiser"],
                    tax_DB = config["tax_DB"]))

    else :
        print("PIRCURST not 'True', will not be in output")


    if config["denoiser"] == "DADA2":
        lst.append(
                expand("{denoiser}/2_denoised/DADA2_denoising_stats.tsv", denoiser = config["denoiser"])
        )


    return lst


## Run the function
rule all:
    input:
        get_final_output(config),



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
