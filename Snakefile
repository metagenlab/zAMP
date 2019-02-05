include: "rules/common_preprocessing/making_sample_dataset.rules"


## A function do define the ouput based on the content of the config file:
def get_final_output(config):

    ## Base output
    lst=[
    #QC
    "QC/multiqc_raw_reads_report.html",
    #"QC/multiqc_DADA2_filtered_reads_report.html",

    ###General output
    expand("{denoiser}/2_denoised/dna-sequences.fasta", denoiser = config["denoiser"]),
    expand("{denoiser}/2_denoised/count-table.qzv", denoiser = config["denoiser"]),
    expand("{denoiser}/2_denoised/rep-seqs.qzv", denoiser = config["denoiser"]),
    expand("{denoiser}/3_classified/rdp/{tax_DB}/dna-sequences_tax_assignments.qzv", denoiser = config["denoiser"], tax_DB = config["tax_DB"]),

    ###notrarefied
    expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/reads/reads_plot_with_filtered.png", denoiser = config["denoiser"], tax_DB = config["tax_DB"]),
    expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/KRONA/{grouping_column}.html", denoiser = config["denoiser"], tax_DB = config["tax_DB"], grouping_column=list(set(all_samples[config["grouping_column"]]))),
    expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/alpha_diversity/{grouping_column}_alpha_divesity.png", denoiser = config["denoiser"], tax_DB = config["tax_DB"], grouping_column=list(set(all_samples[config["grouping_column"]]))),
    expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/ordination/distance_based/{ordination_method_distance_based}/{grouping_column}_d_{ordination_distance}.png", denoiser = config["denoiser"], tax_DB = config["tax_DB"], grouping_column=list(set(all_samples[config["grouping_column"]])), ordination_method_distance_based = config["ordination_method_distance_based"], ordination_distance = config["ordination_distance"]),
    expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/ordination/unconstrained/{ordination_method_unconstrained}/{grouping_column}.png", denoiser = config["denoiser"], tax_DB = config["tax_DB"], grouping_column=list(set(all_samples[config["grouping_column"]])), ordination_method_unconstrained = config["ordination_method_unconstrained"], ordination_distance = config["ordination_distance"]),
    expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/ordination/constrained/{ordination_method_constrained}/{grouping_column}_f_{ordination_factor}.png", denoiser = config["denoiser"], tax_DB = config["tax_DB"], grouping_column=list(set(all_samples[config["grouping_column"]])), ordination_method_constrained = config["ordination_method_constrained"], ordination_factor = config["ordination_factor"]),


    expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/rarefaction_curve.png", denoiser = config["denoiser"], tax_DB = config["tax_DB"]),
    expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/physeq/no_collapse/{filter_tax_rank}_{filter_lineage}_taxfilt_{filter_column_value}_in_{filter_meta_column}_featuresfilt_melted.tsv", denoiser = config["denoiser"], tax_DB = config["tax_DB"], filter_column_value = config["filter_column_value"], filter_meta_column = config["filter_meta_column"], filter_tax_rank = config["filter_tax_rank"], filter_lineage = config["filter_lineage"]),
    expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/physeq/collap_{collapse_level}/{filter_tax_rank}_{filter_lineage}_taxfilt_{filter_column_value}_in_{filter_meta_column}_featuresfilt_melted.tsv", denoiser = config["denoiser"], tax_DB = config["tax_DB"], filter_column_value = config["filter_column_value"], filter_meta_column = config["filter_meta_column"], filter_tax_rank = config["filter_tax_rank"], filter_lineage = config["filter_lineage"] , collapse_level = config["collapse_level"]),

    ###rarefied
    expand("{denoiser}/5_visualization/rdp/{tax_DB}/rarefaction_{rarefaction_value}/alpha_diversity/{grouping_column}_alpha_divesity.png", denoiser = config["denoiser"], tax_DB = config["tax_DB"], grouping_column=list(set(all_samples[config["grouping_column"]])), rarefaction_value=config["rarefaction_value"]),
    expand("{denoiser}/5_visualization/rdp/{tax_DB}/rarefaction_{rarefaction_value}/rarefaction_curve.png", denoiser = config["denoiser"], tax_DB = config["tax_DB"], rarefaction_value=config["rarefaction_value"]),
    expand("{denoiser}/5_visualization/rdp/{tax_DB}/rarefaction_{rarefaction_value}/KRONA/{grouping_column}.html", denoiser = config["denoiser"], tax_DB = config["tax_DB"], grouping_column=list(set(all_samples[config["grouping_column"]])), rarefaction_value=config["rarefaction_value"]),
    expand("{denoiser}/5_visualization/rdp/{tax_DB}/rarefaction_{rarefaction_value}/physeq/collap_{collapse_level}/{filter_tax_rank}_{filter_lineage}_taxfilt_{filter_column_value}_in_{filter_meta_column}_featuresfilt_melted.tsv",  denoiser = config["denoiser"], tax_DB = config["tax_DB"], rarefaction_value=config["rarefaction_value"], filter_column_value = config["filter_column_value"], filter_meta_column = config["filter_meta_column"], filter_tax_rank = config["filter_tax_rank"], filter_lineage = config["filter_lineage"] , collapse_level = config["collapse_level"]),
    expand("{denoiser}/5_visualization/rdp/{tax_DB}/rarefaction_{rarefaction_value}/ordination/distance_based/{ordination_method_distance_based}/{grouping_column}_d_{ordination_distance}.png", denoiser = config["denoiser"], tax_DB = config["tax_DB"], rarefaction_value=config["rarefaction_value"], grouping_column=list(set(all_samples[config["grouping_column"]])), ordination_method_distance_based = config["ordination_method_distance_based"], ordination_distance = config["ordination_distance"]),
    expand("{denoiser}/5_visualization/rdp/{tax_DB}/rarefaction_{rarefaction_value}/ordination/unconstrained/{ordination_method_unconstrained}/{grouping_column}.png", denoiser = config["denoiser"], tax_DB = config["tax_DB"], rarefaction_value=config["rarefaction_value"], grouping_column=list(set(all_samples[config["grouping_column"]])), ordination_method_unconstrained = config["ordination_method_unconstrained"], ordination_distance = config["ordination_distance"]),
    expand("{denoiser}/5_visualization/rdp/{tax_DB}/rarefaction_{rarefaction_value}/ordination/constrained/{ordination_method_constrained}/{grouping_column}_f_{ordination_factor}.png", denoiser = config["denoiser"], tax_DB = config["tax_DB"], rarefaction_value=config["rarefaction_value"],grouping_column=list(set(all_samples[config["grouping_column"]])), ordination_method_constrained = config["ordination_method_constrained"], ordination_factor = config["ordination_factor"]),


    ###ANCOM
    expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/diff_abundance/no_collapse/ANCOM/{filter_tax_rank}_{filter_lineage}_taxfilt_{filter_column_value}_in_{filter_meta_column}_f_{tested_factor}.qzv",  denoiser = config["denoiser"], tax_DB = config["tax_DB"], tested_factor= config["ANCOM_factors"], filter_tax_rank = config["filter_tax_rank"], filter_lineage = config["filter_lineage"], filter_column_value = config["filter_column_value"], filter_meta_column = config["filter_meta_column"]),
    expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/diff_abundance/collap_{collapse_level}/ANCOM/{filter_tax_rank}_{filter_lineage}_taxfilt_{filter_column_value}_in_{filter_meta_column}_f_{tested_factor}.qzv",  denoiser = config["denoiser"], tax_DB = config["tax_DB"], collapse_level = config["collapse_level"], tested_factor= config["ANCOM_factors"], filter_tax_rank = config["filter_tax_rank"], filter_lineage = config["filter_lineage"], filter_column_value = config["filter_column_value"], filter_meta_column = config["filter_meta_column"]),

    ###Gneiss correlation based
    expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/diff_abundance/collap_{collapse_level}/Gneiss/correlation/{filter_tax_rank}_{filter_lineage}_taxfilt_{filter_column_value}_in_{filter_meta_column}/hier_correlation_regression.qzv", denoiser = config["denoiser"], tax_DB = config["tax_DB"], collapse_level = config["collapse_level"], tested_factor= config["ANCOM_factors"], y_balances = list(range(1, 9)), filter_tax_rank = config["filter_tax_rank"], filter_lineage = config["filter_lineage"], filter_column_value = config["filter_column_value"], filter_meta_column = config["filter_meta_column"]),
    expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/diff_abundance/collap_{collapse_level}/Gneiss/correlation/{filter_tax_rank}_{filter_lineage}_taxfilt_{filter_column_value}_in_{filter_meta_column}/hier_correlation_heatmap_{tested_factor}.qzv", denoiser = config["denoiser"], tax_DB = config["tax_DB"], collapse_level = config["collapse_level"], tested_factor= config["ANCOM_factors"], y_balances = list(range(1, 9)), filter_tax_rank = config["filter_tax_rank"], filter_lineage = config["filter_lineage"], filter_column_value = config["filter_column_value"], filter_meta_column = config["filter_meta_column"]),
    expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/diff_abundance/collap_{collapse_level}/Gneiss/correlation/{filter_tax_rank}_{filter_lineage}_taxfilt_{filter_column_value}_in_{filter_meta_column}/hier_correlation_y_{y_balances}_f_{tested_factor}.qzv", denoiser = config["denoiser"], tax_DB = config["tax_DB"], collapse_level = config["collapse_level"], tested_factor= config["ANCOM_factors"], y_balances = list(range(1, 9)), filter_tax_rank = config["filter_tax_rank"], filter_lineage = config["filter_lineage"], filter_column_value = config["filter_column_value"], filter_meta_column = config["filter_meta_column"]),
    expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/diff_abundance/no_collapse/Gneiss/correlation/{filter_tax_rank}_{filter_lineage}_taxfilt_{filter_column_value}_in_{filter_meta_column}/hier_correlation_regression.qzv", denoiser = config["denoiser"], tax_DB = config["tax_DB"], collapse_level = config["collapse_level"], tested_factor= config["ANCOM_factors"], y_balances = list(range(1, 9)), filter_tax_rank = config["filter_tax_rank"], filter_lineage = config["filter_lineage"], filter_column_value = config["filter_column_value"], filter_meta_column = config["filter_meta_column"]),
    expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/diff_abundance/no_collapse/Gneiss/correlation/{filter_tax_rank}_{filter_lineage}_taxfilt_{filter_column_value}_in_{filter_meta_column}/hier_correlation_heatmap_{tested_factor}.qzv", denoiser = config["denoiser"], tax_DB = config["tax_DB"], collapse_level = config["collapse_level"], tested_factor= config["ANCOM_factors"], y_balances = list(range(1, 9)), filter_tax_rank = config["filter_tax_rank"], filter_lineage = config["filter_lineage"], filter_column_value = config["filter_column_value"], filter_meta_column = config["filter_meta_column"]),
    expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/diff_abundance/no_collapse/Gneiss/correlation/{filter_tax_rank}_{filter_lineage}_taxfilt_{filter_column_value}_in_{filter_meta_column}/hier_correlation_y_{y_balances}_f_{tested_factor}.qzv", denoiser = config["denoiser"], tax_DB = config["tax_DB"], collapse_level = config["collapse_level"], tested_factor= config["ANCOM_factors"], y_balances = list(range(1, 9)), filter_tax_rank = config["filter_tax_rank"], filter_lineage = config["filter_lineage"], filter_column_value = config["filter_column_value"], filter_meta_column = config["filter_meta_column"]),

    ###Gneiss phylogeny based
    expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/diff_abundance/collap_{collapse_level}/Gneiss/phylogeny/{filter_tax_rank}_{filter_lineage}_taxfilt_{filter_column_value}_in_{filter_meta_column}/phyl_phylogenetic_regression.qzv", denoiser = config["denoiser"], tax_DB = config["tax_DB"], collapse_level = config["collapse_level"], tested_factor= config["ANCOM_factors"], y_balances = list(range(1, 9)), filter_tax_rank = config["filter_tax_rank"], filter_lineage = config["filter_lineage"], filter_column_value = config["filter_column_value"], filter_meta_column = config["filter_meta_column"]),
    expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/diff_abundance/collap_{collapse_level}/Gneiss/phylogeny/{filter_tax_rank}_{filter_lineage}_taxfilt_{filter_column_value}_in_{filter_meta_column}/phyl_phylogenetic_heatmap_{tested_factor}.qzv", denoiser = config["denoiser"], tax_DB = config["tax_DB"], collapse_level = config["collapse_level"], tested_factor= config["ANCOM_factors"], y_balances = list(range(1, 9)), filter_tax_rank = config["filter_tax_rank"], filter_lineage = config["filter_lineage"], filter_column_value = config["filter_column_value"], filter_meta_column = config["filter_meta_column"]),
    expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/diff_abundance/collap_{collapse_level}/Gneiss/phylogeny/{filter_tax_rank}_{filter_lineage}_taxfilt_{filter_column_value}_in_{filter_meta_column}/phyl_phylogenetic_y_{y_balances}_f_{tested_factor}.qzv", denoiser = config["denoiser"], tax_DB = config["tax_DB"], collapse_level = config["collapse_level"], tested_factor= config["ANCOM_factors"], y_balances = list(range(1, 9)), filter_tax_rank = config["filter_tax_rank"], filter_lineage = config["filter_lineage"], filter_column_value = config["filter_column_value"], filter_meta_column = config["filter_meta_column"]),
    expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/physeq/collap_{collapse_level}/{filter_tax_rank}_{filter_lineage}_taxfilt_{filter_column_value}_in_{filter_meta_column}_featuresfilt/count-table.qzv",  denoiser = config["denoiser"], tax_DB = config["tax_DB"], collapse_level = config["collapse_level"], tested_factor= config["ANCOM_factors"], y_balances = list(range(1, 9)), filter_tax_rank = config["filter_tax_rank"], filter_lineage = config["filter_lineage"], filter_column_value = config["filter_column_value"], filter_meta_column = config["filter_meta_column"]),
    expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/diff_abundance/no_collapse/Gneiss/phylogeny/{filter_tax_rank}_{filter_lineage}_taxfilt_{filter_column_value}_in_{filter_meta_column}/phyl_phylogenetic_regression.qzv", denoiser = config["denoiser"], tax_DB = config["tax_DB"], collapse_level = config["collapse_level"], tested_factor= config["ANCOM_factors"], y_balances = list(range(1, 9)), filter_tax_rank = config["filter_tax_rank"], filter_lineage = config["filter_lineage"], filter_column_value = config["filter_column_value"], filter_meta_column = config["filter_meta_column"]),
    expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/diff_abundance/no_collapse/Gneiss/phylogeny/{filter_tax_rank}_{filter_lineage}_taxfilt_{filter_column_value}_in_{filter_meta_column}/phyl_phylogenetic_heatmap_{tested_factor}.qzv", denoiser = config["denoiser"], tax_DB = config["tax_DB"], collapse_level = config["collapse_level"], tested_factor= config["ANCOM_factors"], y_balances = list(range(1, 9)), filter_tax_rank = config["filter_tax_rank"], filter_lineage = config["filter_lineage"], filter_column_value = config["filter_column_value"], filter_meta_column = config["filter_meta_column"]),
    expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/diff_abundance/no_collapse/Gneiss/phylogeny/{filter_tax_rank}_{filter_lineage}_taxfilt_{filter_column_value}_in_{filter_meta_column}/phyl_phylogenetic_y_{y_balances}_f_{tested_factor}.qzv", denoiser = config["denoiser"], tax_DB = config["tax_DB"], collapse_level = config["collapse_level"], tested_factor= config["ANCOM_factors"], y_balances = list(range(1, 9)), filter_tax_rank = config["filter_tax_rank"], filter_lineage = config["filter_lineage"], filter_column_value = config["filter_column_value"], filter_meta_column = config["filter_meta_column"]),
    expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/physeq/no_collapse/{filter_tax_rank}_{filter_lineage}_taxfilt_{filter_column_value}_in_{filter_meta_column}_featuresfilt/count-table.qzv",  denoiser = config["denoiser"], tax_DB = config["tax_DB"], collapse_level = config["collapse_level"], tested_factor= config["ANCOM_factors"], y_balances = list(range(1, 9)), filter_tax_rank = config["filter_tax_rank"], filter_lineage = config["filter_lineage"], filter_column_value = config["filter_column_value"], filter_meta_column = config["filter_meta_column"])

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


