include: "rules/common_preprocessing/making_sample_dataset.rules"

rule all:
    input:
        # QC
        "QC/multiqc_raw_reads_report.html",
        "QC/multiqc_DADA2_filtered_reads_report.html",

        # DADA2
        #"DADA2/2_denoised/dna-sequences.fasta",
        #"DADA2/2_denoised/DADA2_denoising_stats.tsv",
        #"DADA2/2_denoised/count-table.qzv",
        #"DADA2/2_denoised/rep-seqs.qzv",
        #"DADA2/5_visualization/rdp/ezbiocloud_valentin/phyloseq_melted_table.tsv",
        expand("DADA2/5_visualization/rdp/ezbiocloud_valentin/KRONA/{grouping_column}.html", grouping_column=list(set(all_samples[config["grouping_column"]]))),
        expand("DADA2/5_visualization/rdp/ezbiocloud_marta/KRONA/{grouping_column}.html", grouping_column=list(set(all_samples[config["grouping_column"]]))),

        # vsearch
        #"vsearch/2_denoised/all_samples_reads_count.txt",
        #"vsearch/2_denoised/count-table.qzv",
        #"vsearch/2_denoised/rep-seqs.qzv",
        #"vsearch/3_classified/rdp/ezbiocloud_valentin/otus_tax_table.txt",
        #"vsearch/5_visualization/rdp/ezbiocloud_valentin/phyloseq_melted_table.tsv",
        expand("vsearch/5_visualization/rdp/ezbiocloud_marta/KRONA/{grouping_column}.html", grouping_column=list(set(all_samples[config["grouping_column"]]))),
        expand("vsearch/5_visualization/rdp/ezbiocloud_valentin/KRONA/{grouping_column}.html", grouping_column=list(set(all_samples[config["grouping_column"]]))),

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
include: "rules/common_visualization/General_plotting.rules"

