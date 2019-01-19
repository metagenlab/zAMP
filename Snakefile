include: "rules/common_preprocessing/making_sample_dataset.rules"

rule all:
    input:
        #QC
        "QC/multiqc_raw_reads_report.html",
        #"QC/multiqc_DADA2_filtered_reads_report.html",

        #DADA2
        "DADA2/2_denoised/dna-sequences.fasta",
        "DADA2/2_denoised/count-table.qzv",
        "DADA2/2_denoised/rep-seqs.qzv",
        "DADA2/2_denoised/DADA2_denoising_stats.tsv",
        #rarefied
        expand("DADA2/5_visualization/rdp/ezbiocloud_valentin/rarefaction_{rarefaction_value}/alpha_diversity/{grouping_column}_alpha_divesity.png", grouping_column=list(set(all_samples[config["grouping_column"]])), rarefaction_value=config["rarefaction_value"]),
        expand("DADA2/5_visualization/rdp/ezbiocloud_valentin/rarefaction_{rarefaction_value}/ordination/{grouping_column}_{ordination_distance}.png",grouping_column=list(set(all_samples[config["grouping_column"]])), rarefaction_value=config["rarefaction_value"], ordination_distance = config["ordination_distance"]),
        expand("DADA2/5_visualization/rdp/ezbiocloud_valentin/rarefaction_{rarefaction_value}/rarefaction_curve.png", rarefaction_value=config["rarefaction_value"]),
        expand("DADA2/5_visualization/rdp/ezbiocloud_valentin/rarefaction_{rarefaction_value}/phyloseq_object", rarefaction_value=config["rarefaction_value"]),
        expand("DADA2/5_visualization/rdp/ezbiocloud_valentin/rarefaction_{rarefaction_value}/KRONA/{grouping_column}.html", grouping_column=list(set(all_samples[config["grouping_column"]])), rarefaction_value=config["rarefaction_value"]),
        #notrarefied
        "DADA2/5_visualization/rdp/ezbiocloud_valentin/norarefaction/rarefaction_curve.png",
        "DADA2/5_visualization/rdp/ezbiocloud_valentin/norarefaction/reads/reads_plot_with_filtered.png",
        expand("DADA2/5_visualization/rdp/ezbiocloud_marta/norarefaction/KRONA/{grouping_column}.html", grouping_column=list(set(all_samples[config["grouping_column"]]))),
        expand("DADA2/5_visualization/rdp/ezbiocloud_valentin/norarefaction/KRONA/{grouping_column}.html", grouping_column=list(set(all_samples[config["grouping_column"]]))),
        expand("DADA2/5_visualization/rdp/ezbiocloud_valentin/norarefaction/alpha_diversity/{grouping_column}_alpha_divesity.png", grouping_column=list(set(all_samples[config["grouping_column"]]))),
        expand("DADA2/5_visualization/rdp/ezbiocloud_valentin/norarefaction/ordination/{grouping_column}_{ordination_distance}.png", grouping_column=list(set(all_samples[config["grouping_column"]])), ordination_distance = config["ordination_distance"]),
        "DADA2/5_visualization/rdp/ezbiocloud_valentin/norarefaction/rarefaction_curve.png",
        "DADA2/5_visualization/rdp/ezbiocloud_valentin/norarefaction/phyloseq_object",

        #expand("DADA2/5_visualization/rdp/ezbiocloud_valentin/norarefaction/adundancy_comparison/taxa_collapse_{collapse_level}/ANCOM_{ANCOM_factor}.qzv", collapse_level = config["collapse_level"], ANCOM_factor= config["ANCOM_factor"]),
        #"DADA2/3_classified/rdp/ezbiocloud_valentin/dna-sequences_tax_assignments.qzv",
        #"DADA2/5_visualization/rdp/ezbiocloud_valentin/norarefaction/adundancy_comparison/filtered_samples.qzv",

        expand("DADA2/5_visualization/rdp/ezbiocloud_valentin/norarefaction/ANCOM/taxa_collapse_{collapse_level}/ANCOM_{ANCOM_factor}.qzv", collapse_level = config["collapse_level"], ANCOM_factor= config["ANCOM_factor"]),
        "DADA2/3_classified/rdp/ezbiocloud_valentin/dna-sequences_tax_assignments.qzv",
        "DADA2/5_visualization/rdp/ezbiocloud_valentin/norarefaction/ANCOM/filtered_samples.qzv",

        #vsearch
        "vsearch/2_denoised/all_samples_reads_count.txt",
        "vsearch/2_denoised/count-table.qzv",
        "vsearch/2_denoised/rep-seqs.qzv",
        "vsearch/3_classified/rdp/ezbiocloud_valentin/otus_tax_table.txt",
        #rarefied
        expand("vsearch/5_visualization/rdp/ezbiocloud_valentin/rarefaction_{rarefaction_value}/alpha_diversity/{grouping_column}_alpha_divesity.png", grouping_column=list(set(all_samples[config["grouping_column"]])), rarefaction_value=config["rarefaction_value"]),
        expand("vsearch/5_visualization/rdp/ezbiocloud_valentin/rarefaction_{rarefaction_value}/ordination/{grouping_column}_{ordination_distance}.png",grouping_column=list(set(all_samples[config["grouping_column"]])), rarefaction_value=config["rarefaction_value"], ordination_distance = config["ordination_distance"]),
        expand("vsearch/5_visualization/rdp/ezbiocloud_valentin/rarefaction_{rarefaction_value}/rarefaction_curve.png", rarefaction_value=config["rarefaction_value"]),
        expand("vsearch/5_visualization/rdp/ezbiocloud_valentin/rarefaction_{rarefaction_value}/phyloseq_object", rarefaction_value=config["rarefaction_value"]),
        expand("vsearch/5_visualization/rdp/ezbiocloud_valentin/rarefaction_{rarefaction_value}/KRONA/{grouping_column}.html", grouping_column=list(set(all_samples[config["grouping_column"]])), rarefaction_value=config["rarefaction_value"]),
        #notrarefied
        #"vsearch/5_visualization/rdp/ezbiocloud_valentin/norarefaction/rarefaction_curve.png",
        "vsearch/5_visualization/rdp/ezbiocloud_valentin/norarefaction/reads/reads_plot_with_filtered.png",
        expand("vsearch/5_visualization/rdp/ezbiocloud_marta/norarefaction/KRONA/{grouping_column}.html", grouping_column=list(set(all_samples[config["grouping_column"]]))),
        expand("vsearch/5_visualization/rdp/ezbiocloud_valentin/norarefaction/KRONA/{grouping_column}.html", grouping_column=list(set(all_samples[config["grouping_column"]]))),
        expand("vsearch/5_visualization/rdp/ezbiocloud_valentin/norarefaction/alpha_diversity/{grouping_column}_alpha_divesity.png", grouping_column=list(set(all_samples[config["grouping_column"]]))),
        expand("vsearch/5_visualization/rdp/ezbiocloud_valentin/norarefaction/ordination/{grouping_column}_bray.png",grouping_column=list(set(all_samples[config["grouping_column"]])), ordination_distance = config["ordination_distance"]),
        "vsearch/5_visualization/rdp/ezbiocloud_valentin/norarefaction/rarefaction_curve.png",
        "vsearch/5_visualization/rdp/ezbiocloud_valentin/norarefaction/phyloseq_object",


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
include: "rules/common_visualization/to_Phyloseq.rules"
include: "rules/common_visualization/rarefy.rules"
include: "rules/common_visualization/General_plotting.rules"
include: "rules/Qiime2_stat_analysis/Qiime2_stat_analysis.rules"
