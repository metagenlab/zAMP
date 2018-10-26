rule all:
    input:
         "DADA2/2_denoised/dna-sequences.fasta",
         "DADA2/2_denoised/DADA2_denoising_stats.tsv",
        #"qiime2/2_denoised/deblur/seqs.qzv",
        #"qiime2/2_denoised/deblur/table.qzv",
        #"qiime2/2_denoised/deblur/stats.qzv",
        #"qiime2/1a_cutadapt_primers_trimming/trimmed_sequences.qzv",
        #"qiime2/1c_q_score_filtering/quality_filtered_stats.qzv",
        #"qiime2/3_classified/deblur/RDP/ezbiocloud_valentin/dna-sequences_tax_assignments.txt",
        #"qiime2/4_tree/rooted-tree.qza",
         "DADA2/3_classified/rdp/ezbiocloud_marta/dna-sequences_tax_assignments.txt",
         "DADA2/4_tree/rooted-tree.qza",
         "DADA2/2_denoised/rep-seqs.qzv",
         "DADA2/2_denoised/count-table.qzv",
         "QC/multiqc_filtered_reads_report.html",
         "QC/multiqc_raw_reads_report.html",
         "DADA2/5_visualization/rdp/ezbiocloud_marta/phyloseq_object",
         "DADA2/5_visualization/rdp/ezbiocloud_marta/phyloseq_melted_table.tsv",
         expand("{tool}/5_visualization/{classifier}/{db_taxonomy}/KRONA/{sample}", sample=list(read_naming.keys()))



include: "rules/get_reads/making_sample_dataset.rules"
include: "rules/get_reads/get_reads.rules"
include: "rules/get_reads/get_sras.rules"
# include: "rules/0_importing_data.rules"
include: "rules/QC/QC.rules"
include: "rules/DADA2/1_trim_for_DADA2.rules"
include: "rules/DADA2/2_DADA2.rules"
include: "rules/QIIME/3_assign_taxonomy.rules"
include: "rules/QIIME2/4_tree.rules"
include: "rules/QIIME2/Import_to_QIIME2.rules"
# include: "rules/qiime.rules"
include: "rules/Visualization/DADA2_to_Phyloseq.rules"
include: "rules/Visualization/General_plotting.rules"
