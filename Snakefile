include: "rules/common_preprocessing/making_sample_dataset.rules"

rule all:
    input:
        #QC
        "QC/multiqc_raw_reads_report.html",
        #"QC/multiqc_DADA2_filtered_reads_report.html",

        #DADA2
        "DADA2/2_denoised/count-table.qzv",
        "DADA2/2_denoised/rep-seqs.qzv",
        #rarefied
        expand("DADA2/5_visualization/rdp/ezbiocloud_valentin/rarefaction_" + str(config["rarefaction_value"]) + "/alpha_diversity/{grouping_column}_alpha_divesity.png", grouping_column=list(set(all_samples[config["grouping_column"]]))),
        "DADA2/5_visualization/rdp/ezbiocloud_valentin/rarefaction_" + str(config["rarefaction_value"]) + "/ordination/unifrac.png",
        "DADA2/5_visualization/rdp/ezbiocloud_valentin/rarefaction_" + str(config["rarefaction_value"]) + "/rarefaction_curve.png",
        "DADA2/5_visualization/rdp/ezbiocloud_valentin/rarefaction_" + str(config["rarefaction_value"]) + "/phyloseq_object",
        #notrarefied
        "DADA2/5_visualization/rdp/ezbiocloud_valentin/norarefaction/rarefaction_curve.png",
        "DADA2/5_visualization/rdp/ezbiocloud_valentin/norarefaction/reads/reads_plot_with_filtered.png",
        #expand("DADA2/5_visualization/rdp/ezbiocloud_marta/norarefaction/KRONA/{grouping_column}.html", grouping_column=list(set(all_samples[config["grouping_column"]]))),
        expand("DADA2/5_visualization/rdp/ezbiocloud_valentin/norarefaction/KRONA/{grouping_column}.html", grouping_column=list(set(all_samples[config["grouping_column"]]))),
        expand("DADA2/5_visualization/rdp/ezbiocloud_valentin/norarefaction/alpha_diversity/{grouping_column}_alpha_divesity.png", grouping_column=list(set(all_samples[config["grouping_column"]]))),
        expand("DADA2/5_visualization/rdp/ezbiocloud_valentin/norarefaction/ordination/{grouping_column}_unifrac.png", grouping_column=list(set(all_samples[config["grouping_column"]]))),
        "DADA2/5_visualization/rdp/ezbiocloud_valentin/norarefaction/rarefaction_curve.png",
        "DADA2/5_visualization/rdp/ezbiocloud_valentin/norarefaction/phyloseq_object",



        #vsearch
        "vsearch/2_denoised/all_samples_reads_count.txt",
        "vsearch/2_denoised/count-table.qzv",
        "vsearch/2_denoised/rep-seqs.qzv",
        "vsearch/3_classified/rdp/ezbiocloud_valentin/otus_tax_table.txt",
        #rarefied
        expand("vsearch/5_visualization/rdp/ezbiocloud_valentin/rarefaction_" + str(config["rarefaction_value"]) + "/alpha_diversity/{grouping_column}_alpha_divesity.png", grouping_column=list(set(all_samples[config["grouping_column"]]))),
        "vsearch/5_visualization/rdp/ezbiocloud_valentin/rarefaction_" + str(config["rarefaction_value"]) + "/ordination/unifrac.png",
        "vsearch/5_visualization/rdp/ezbiocloud_valentin/rarefaction_" + str(config["rarefaction_value"]) + "/rarefaction_curve.png",
        "vsearch/5_visualization/rdp/ezbiocloud_valentin/rarefaction_" + str(config["rarefaction_value"]) + "/phyloseq_object",
        #notrarefied
        "vsearch/5_visualization/rdp/ezbiocloud_valentin/norarefaction/rarefaction_curve.png",
        "vsearch/5_visualization/rdp/ezbiocloud_valentin/norarefaction/reads/reads_plot_with_filtered.png",
        #expand("vsearch/5_visualization/rdp/ezbiocloud_marta/norarefaction/KRONA/{grouping_column}.html", grouping_column=list(set(all_samples[config["grouping_column"]]))),
        expand("vsearch/5_visualization/rdp/ezbiocloud_valentin/norarefaction/KRONA/{grouping_column}.html", grouping_column=list(set(all_samples[config["grouping_column"]]))),
        expand("vsearch/5_visualization/rdp/ezbiocloud_valentin/norarefaction/alpha_diversity/{grouping_column}_alpha_divesity.png", grouping_column=list(set(all_samples[config["grouping_column"]]))),
        expand("vsearch/5_visualization/rdp/ezbiocloud_valentin/norarefaction/ordination/{grouping_column}_unifrac.png",grouping_column=list(set(all_samples[config["grouping_column"]]))),
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
