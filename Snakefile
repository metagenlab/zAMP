include: "rules/common_preprocessing/making_sample_dataset.rules"

rule all:
    input:
        #QC
        "QC/multiqc_raw_reads_report.html",
        #"QC/multiqc_DADA2_filtered_reads_report.html",

        #DADA2
        expand("{denoiser}/2_denoised/dna-sequences.fasta", denoiser = config["denoiser"]),
        expand("{denoiser}/2_denoised/count-table.qzv", denoiser = config["denoiser"]),
        expand("{denoiser}/2_denoised/rep-seqs.qzv", denoiser = config["denoiser"]),
        #expand("{denoiser}/2_denoised/DADA2_denoising_stats.tsv", denoiser = config["denoiser"]),
        expand("{denoiser}/3_classified/rdp/{tax_DB}/dna-sequences_tax_assignments.qzv", denoiser = config["denoiser"], tax_DB = config["tax_DB"]),
        #rarefied
        expand("{denoiser}/5_visualization/rdp/{tax_DB}/rarefaction_{rarefaction_value}/alpha_diversity/{grouping_column}_alpha_divesity.png", denoiser = config["denoiser"], tax_DB = config["tax_DB"], grouping_column=list(set(all_samples[config["grouping_column"]])), rarefaction_value=config["rarefaction_value"]),
        expand("{denoiser}/5_visualization/rdp/{tax_DB}/rarefaction_{rarefaction_value}/ordination/{grouping_column}_{ordination_distance}.png", denoiser = config["denoiser"], tax_DB = config["tax_DB"],grouping_column=list(set(all_samples[config["grouping_column"]])), rarefaction_value=config["rarefaction_value"], ordination_distance = config["ordination_distance"]),
        #expand("{denoiser}/5_visualization/rdp/{tax_DB}/rarefaction_{rarefaction_value}/rarefaction_curve.png", denoiser = config["denoiser"], tax_DB = config["tax_DB"], rarefaction_value=config["rarefaction_value"]),
        #expand("{denoiser}/5_visualization/rdp/{tax_DB}/rarefaction_{rarefaction_value}/phyloseq_object", denoiser = config["denoiser"], tax_DB = config["tax_DB"], rarefaction_value=config["rarefaction_value"]),
        expand("{denoiser}/5_visualization/rdp/{tax_DB}/rarefaction_{rarefaction_value}/KRONA/{grouping_column}.html", denoiser = config["denoiser"], tax_DB = config["tax_DB"], grouping_column=list(set(all_samples[config["grouping_column"]])), rarefaction_value=config["rarefaction_value"]),
        #notrarefied
        expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/rarefaction_curve.png", denoiser = config["denoiser"], tax_DB = config["tax_DB"]),
        expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/reads/reads_plot_with_filtered.png", denoiser = config["denoiser"], tax_DB = config["tax_DB"]),
        expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/KRONA/{grouping_column}.html", denoiser = config["denoiser"], tax_DB = config["tax_DB"], grouping_column=list(set(all_samples[config["grouping_column"]]))),
        expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/KRONA/{grouping_column}.html", denoiser = config["denoiser"], tax_DB = config["tax_DB"], grouping_column=list(set(all_samples[config["grouping_column"]]))),
        expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/alpha_diversity/{grouping_column}_alpha_divesity.png", denoiser = config["denoiser"], tax_DB = config["tax_DB"], grouping_column=list(set(all_samples[config["grouping_column"]]))),
        expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/ordination/{grouping_column}_{ordination_distance}.png", denoiser = config["denoiser"], tax_DB = config["tax_DB"], grouping_column=list(set(all_samples[config["grouping_column"]])), ordination_distance = config["ordination_distance"]),
        expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/rarefaction_curve.png", denoiser = config["denoiser"], tax_DB = config["tax_DB"]),
        expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/physeq/{filter_tax_rank}_{filter_lineage}_taxfilt_{column_value}_in_{meta_column}_featuresfilt_{collapse_level}_collapsed_melted.tsv", denoiser = config["denoiser"], tax_DB = config["tax_DB"], column_value = config["column_value"], meta_column = config["meta_column"], filter_tax_rank = config["filter_tax_rank"], filter_lineage = config["filter_lineage"] , collapse_level = config["collapse_level"]),


        #ANCOM
        #expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/differential_abundance/no_collapse/ANCOM/ANCOM_{tested_factor}.qzv", denoiser = config["denoiser"], tax_DB = config["tax_DB"], tested_factor= config["ANCOM_factors"]),
        #expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/differential_abundance/taxa_collapse_{collapse_level}/ANCOM/ANCOM_{tested_factor}.qzv", denoiser = config["denoiser"], tax_DB = config["tax_DB"], collapse_level = config["collapse_level"], tested_factor= config["ANCOM_factors"]),
        #expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/differential_abundance/no_collapse/filtered_samples.qzv", denoiser = config["denoiser"], tax_DB = config["tax_DB"]),

        #Gneiss
        #expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/differential_abundance/no_collapse/Gneiss/correlation/regression_hier_correlation.qzv", denoiser = config["denoiser"], tax_DB = config["tax_DB"]),
        #expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/differential_abundance/no_collapse/Gneiss/gradient/regression_hier_{tested_factor}.qzv", denoiser = config["denoiser"], tax_DB = config["tax_DB"], tested_factor = config["Gneiss_gradient_clustering"]),
        #expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/differential_abundance/taxa_collapse_{collapse_level}/Gneiss/correlation/regression_hier_correlation.qzv", denoiser = config["denoiser"], tax_DB = config["tax_DB"], collapse_level = config["collapse_level"]),
        #expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/differential_abundance/taxa_collapse_{collapse_level}/Gneiss/gradient/regression_hier_{clustering_factor}.qzv", denoiser = config["denoiser"], tax_DB = config["tax_DB"], clustering_factor = config["Gneiss_gradient_clustering"],collapse_level = config["collapse_level"]),
        #expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/differential_abundance/taxa_collapse_{collapse_level}/Gneiss/correlation/{tested_factor}_heatmap_hier_correlation.qzv", denoiser = config["denoiser"], tax_DB = config["tax_DB"], collapse_level = config["collapse_level"], tested_factor= config["ANCOM_factors"]),
        #expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/differential_abundance/no_collapse/Gneiss/correlation/{tested_factor}_heatmap_hier_correlation.qzv", denoiser = config["denoiser"], tax_DB = config["tax_DB"], tested_factor= config["ANCOM_factors"]),
        #expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/differential_abundance/taxa_collapse_{collapse_level}/Gneiss/correlation/{tested_factor}_taxbalance_hier_correlation_y_{y_balances}.qzv", denoiser = config["denoiser"], tax_DB = config["tax_DB"], collapse_level = config["collapse_level"], tested_factor= config["ANCOM_factors"], y_balances = list(range(1, 9))),
        #expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/differential_abundance/no_collapse/Gneiss/phylogeny/regression_phyl_phylogenetic.qzv", denoiser = config["denoiser"], tax_DB = config["tax_DB"]),
        #expand("{denoiser}/5_visualization/rdp/{tax_DB}/norarefaction/differential_abundance/taxa_collapse_{collapse_level}/Gneiss/phylogeny/regression_phyl_phylogenetic.qzv", denoiser = config["denoiser"], tax_DB = config["tax_DB"], collapse_level = config["collapse_level"]),



        #vsearch
        #"vsearch/2_denoised/all_samples_reads_count.txt",
        #"vsearch/2_denoised/count-table.qzv",
        #"vsearch/2_denoised/rep-seqs.qzv",
        #"vsearch/3_classified/rdp/ezbiocloud_valentin/otus_tax_table.txt",
        #rarefied
        #expand("vsearch/5_visualization/rdp/ezbiocloud_valentin/rarefaction_{rarefaction_value}/alpha_diversity/{grouping_column}_alpha_divesity.png", grouping_column=list(set(all_samples[config["grouping_column"]])), rarefaction_value=config["rarefaction_value"]),
        #expand("vsearch/5_visualization/rdp/ezbiocloud_valentin/rarefaction_{rarefaction_value}/ordination/{grouping_column}_{ordination_distance}.png",grouping_column=list(set(all_samples[config["grouping_column"]])), rarefaction_value=config["rarefaction_value"], ordination_distance = config["ordination_distance"]),
        #expand("vsearch/5_visualization/rdp/ezbiocloud_valentin/rarefaction_{rarefaction_value}/rarefaction_curve.png", rarefaction_value=config["rarefaction_value"]),
        #expand("vsearch/5_visualization/rdp/ezbiocloud_valentin/rarefaction_{rarefaction_value}/phyloseq_object", rarefaction_value=config["rarefaction_value"]),
        #expand("vsearch/5_visualization/rdp/ezbiocloud_valentin/rarefaction_{rarefaction_value}/KRONA/{grouping_column}.html", grouping_column=list(set(all_samples[config["grouping_column"]])), rarefaction_value=config["rarefaction_value"]),
        #notrarefied
        #"vsearch/5_visualization/rdp/ezbiocloud_valentin/norarefaction/rarefaction_curve.png",
        #"vsearch/5_visualization/rdp/ezbiocloud_valentin/norarefaction/reads/reads_plot_with_filtered.png",
        #expand("vsearch/5_visualization/rdp/ezbiocloud_marta/norarefaction/KRONA/{grouping_column}.html", grouping_column=list(set(all_samples[config["grouping_column"]]))),
        #expand("vsearch/5_visualization/rdp/ezbiocloud_valentin/norarefaction/KRONA/{grouping_column}.html", grouping_column=list(set(all_samples[config["grouping_column"]]))),
        #expand("vsearch/5_visualization/rdp/ezbiocloud_valentin/norarefaction/alpha_diversity/{grouping_column}_alpha_divesity.png", grouping_column=list(set(all_samples[config["grouping_column"]]))),
        #expand("vsearch/5_visualization/rdp/ezbiocloud_valentin/norarefaction/ordination/{grouping_column}_bray.png",grouping_column=list(set(all_samples[config["grouping_column"]])), ordination_distance = config["ordination_distance"]),
        #"vsearch/5_visualization/rdp/ezbiocloud_valentin/norarefaction/rarefaction_curve.png",
        #"vsearch/5_visualization/rdp/ezbiocloud_valentin/norarefaction/phyloseq_object",


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
include: "rules/common_visualization/rarefy.rules"
include: "rules/common_visualization/General_plotting.rules"
include: "rules/Qiime2_stat_analysis/Qiime2_stat_analysis.rules"
