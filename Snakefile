include: "rules/common_preprocessing/making_sample_dataset.rules"


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
    file_list = []

    for i in set(column_of_interest):
        combined_values = expand("{column_of_interest}/{column_values}", column_of_interest = i, column_values = list(set(all_samples[i])))
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

    ## Rarefied visualiation

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
        expand("{denoiser}/5_visualization/rdp/{tax_DB}/{rarefaction_value}/physeq/{collapse_key}/2_filter_samples/{filter_tax_rank}_{filter_lineage}_{filter_column_value}_in_{filter_meta_column}_trsf_melted.tsv",
            denoiser = config["denoiser"],
            tax_DB = config["tax_DB"],
            collapse_key = get_taxa_collapse_level_key(config["collapse_level"]),
            rarefaction_value = get_rarefaction_key(config["rarefaction_value"]),
            filter_tax_rank = config["filter_tax_rank"],
            filter_lineage = config["filter_lineage"],
            filter_column_value = config["filter_column_value"],
            filter_meta_column = config["filter_meta_column"]),


    ## Statistical analyses

        ### ANCOM
            expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/diff_abundance/{collapse_key}/ANCOM/{filter_tax_rank}_{filter_lineage}_taxfilt_{filter_column_value}_in_{filter_meta_column}_f_{tested_factor}.qzv",
                denoiser = config["denoiser"],
                tax_DB = config["tax_DB"],
                collapse_key = get_taxa_collapse_level_key(config["collapse_level"]),
                filter_tax_rank = config["filter_tax_rank"],
                filter_lineage = config["filter_lineage"],
                filter_column_value = config["filter_column_value"],
                filter_meta_column = config["filter_meta_column"],
                tested_factor= config["ANCOM_factors"]),


        ### Gneiss
            #### Correlation based
                ##### Regression - Taxa collapse
            expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/diff_abundance/{collapse_key}/Gneiss/correlation/{filter_tax_rank}_{filter_lineage}_taxfilt_{filter_column_value}_in_{filter_meta_column}/hier_correlation_regression.qzv",
                denoiser = config["denoiser"],
                tax_DB = config["tax_DB"],
                collapse_key = get_taxa_collapse_level_key(config["collapse_level"]),
                filter_tax_rank = config["filter_tax_rank"],
                filter_lineage = config["filter_lineage"],
                filter_column_value = config["filter_column_value"],
                filter_meta_column = config["filter_meta_column"]),


                ##### Heatmaps - Taxa collapse
            expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/diff_abundance/{collapse_key}/Gneiss/correlation/{filter_tax_rank}_{filter_lineage}_taxfilt_{filter_column_value}_in_{filter_meta_column}/hier_correlation_heatmap_{tested_factor}.qzv",
                denoiser = config["denoiser"],
                tax_DB = config["tax_DB"],
                collapse_key = get_taxa_collapse_level_key(config["collapse_level"]),
                filter_tax_rank = config["filter_tax_rank"],
                filter_lineage = config["filter_lineage"],
                filter_column_value = config["filter_column_value"],
                filter_meta_column = config["filter_meta_column"],
                tested_factor= config["ANCOM_factors"]),


                ##### Balances - Taxa collapase
            expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/diff_abundance/{collapse_key}/Gneiss/correlation/{filter_tax_rank}_{filter_lineage}_taxfilt_{filter_column_value}_in_{filter_meta_column}/hier_correlation_y_{y_balances}_f_{tested_factor}.qzv",
                denoiser = config["denoiser"],
                tax_DB = config["tax_DB"],
                collapse_key = get_taxa_collapse_level_key(config["collapse_level"]),
                filter_tax_rank = config["filter_tax_rank"],
                filter_lineage = config["filter_lineage"],
                filter_column_value = config["filter_column_value"],
                filter_meta_column = config["filter_meta_column"],
                tested_factor= config["ANCOM_factors"],
                y_balances = list(range(1, 9))),


            ### Phylogeny based
                ##### Regression - Taxa collapse
            expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/diff_abundance/{collapse_key}/Gneiss/phylogeny/{filter_tax_rank}_{filter_lineage}_taxfilt_{filter_column_value}_in_{filter_meta_column}/phyl_phylogenetic_regression.qzv",
                denoiser = config["denoiser"],
                tax_DB = config["tax_DB"],
                collapse_key = get_taxa_collapse_level_key(config["collapse_level"]),
                filter_tax_rank = config["filter_tax_rank"],
                filter_lineage = config["filter_lineage"],
                filter_column_value = config["filter_column_value"],
                filter_meta_column = config["filter_meta_column"]),


                ##### Heatmaps - Taxa collapse
            expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/diff_abundance/{collapse_key}/Gneiss/phylogeny/{filter_tax_rank}_{filter_lineage}_taxfilt_{filter_column_value}_in_{filter_meta_column}/phyl_phylogenetic_heatmap_{tested_factor}.qzv",
                denoiser = config["denoiser"],
                tax_DB = config["tax_DB"],
                collapse_key = get_taxa_collapse_level_key(config["collapse_level"]),
                y_balances = list(range(1, 9)),
                filter_tax_rank = config["filter_tax_rank"],
                filter_lineage = config["filter_lineage"],
                filter_column_value = config["filter_column_value"],
                filter_meta_column = config["filter_meta_column"],
                tested_factor = config["ANCOM_factors"]),


                ##### Balances - Taxa collapase
            expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/diff_abundance/{collapse_key}/Gneiss/phylogeny/{filter_tax_rank}_{filter_lineage}_taxfilt_{filter_column_value}_in_{filter_meta_column}/phyl_phylogenetic_y_{y_balances}_f_{tested_factor}.qzv",
                denoiser = config["denoiser"],
                tax_DB = config["tax_DB"],
                collapse_key = get_taxa_collapse_level_key(config["collapse_level"]),
                filter_tax_rank = config["filter_tax_rank"],
                filter_lineage = config["filter_lineage"],
                filter_column_value = config["filter_column_value"],
                filter_meta_column = config["filter_meta_column"],
                tested_factor= config["ANCOM_factors"],
                y_balances = list(range(1, 9))),

    ###Gneiss gradient based to be implemented
    ]

    ## Conditional output
    if config["denoiser"] == "DADA2":
        lst.append(
                expand("{denoiser}/2_denoised/DADA2_denoising_stats.tsv", denoiser = config["denoiser"]))
    return lst


## Run the function
rule all:
    input:
        #QC
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
include: "rules/common_visualization/Qiime2_stat_analysis_phyloseq_preprocessing.rules"


