include: "rules/get_reads/making_sample_dataset.rules"

rule all:
    input:
        "DADA2/2_denoised/dna-sequences.fasta",
        "DADA2/2_denoised/DADA2_denoising_stats.tsv",
        "DADA2/3_classified/rdp/ezbiocloud_marta/dna-sequences_tax_assignments.txt",
        "DADA2/4_tree/rooted-tree.qza",
        "DADA2/2_denoised/rep-seqs.qzv",
        "DADA2/2_denoised/count-table.qzv",
        "QC/multiqc_filtered_reads_report.html",
        "QC/multiqc_raw_reads_report.html",
        "DADA2/5_visualization/rdp/ezbiocloud_marta/phyloseq_object",
        "DADA2/5_visualization/rdp/ezbiocloud_marta/phyloseq_melted_table.tsv",
        expand("DADA2/5_visualization/rdp/ezbiocloud_marta/KRONA/{sample}.html", sample=list(read_naming.keys())),
        "Marta/3_classified/rdp/ezbiocloud_marta/dna-sequences_tax_assignments.txt",
        "Marta/3_classified/rdp/ezbiocloud_marta/dna-sequences_tax_assignments_reformatted.txt",
        "Marta/2_denoised/count_table.txt",
        "Marta/2_denoised/all_samples_reads_count.txt"

include: "rules/get_reads/get_reads.rules"
include: "rules/get_reads/get_sras.rules"
include: "rules/QC/QC.rules"
include: "rules/DADA2/1_trim_for_DADA2.rules"
include: "rules/DADA2/2_DADA2.rules"
include: "rules/QIIME/3_assign_taxonomy.rules"
include: "rules/QIIME2/4_tree.rules"
include: "rules/QIIME2/Import_to_QIIME2.rules"
include: "rules/Visualization/DADA2_to_Phyloseq.rules"
include: "rules/Visualization/General_plotting.rules"
include: "rules/Martapipeline/1_reads_processing.rules"

