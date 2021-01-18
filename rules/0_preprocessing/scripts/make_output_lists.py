## Define lits of outputs as that will be called by function from the Snakefile depending from the values in config.

### Import list of functions to handle the output
include: "make_output_fcts.py"

### Basic output
#### MuliQC
MultiQC = expand("QC/{RUN}_multiqc_raw_reads_report.html", RUN = set(all_samples[config["run_column"]]))
MultiQC.append("QC/multiqc_raw_reads_report.html")


#### Light pipeline output
light_output = expand("{denoiser}/2_denoised/dna-sequences.fasta", denoiser = config["denoiser"])
light_output.append(expand("{denoiser}/2_denoised/count_table.tsv", denoiser = config["denoiser"]))
light_output.append(expand("{denoiser}/3_classified/{classifier}_{tax_DB}/dna-sequences_tax_assignments.qzv", classifier = config["classifier"], denoiser = config["denoiser"], tax_DB = config["tax_DB_name"]))

#### DADA2 stat table
DADA2_stats_tables = "DADA2/2_denoised/DADA2_denoising_stats.tsv"


### Basic evaluation plots
basic_plots = expand("{denoiser}/5_visualization/{classifier}_{tax_DB}/reads/reads_plot_with_filtered.png",
                                 classifier = config["classifier"],
                                 denoiser = config["denoiser"],
                                 tax_DB = config["tax_DB_name"])



basic_plots.append(expand("{denoiser}/5_visualization/{classifier}_{tax_DB}/reads/reads_plot_with_{filter_lineage}_in_{filter_tax_rank}_filtered.png",
                                denoiser = config["denoiser"],
                                classifier = config["classifier"],
                                tax_DB = config["tax_DB_name"],
                                filter_tax_rank = config["filter_tax_rank"],
                                filter_lineage = config["filter_lineage"]))


basic_plots.append(expand("{denoiser}/5_visualization/{classifier}_{tax_DB}/rarefaction_curve/nofiltering/{rarefaction_value}_rarefaction_curve_{grouping_column}.png",
                                 denoiser = config["denoiser"],
                                 classifier = config["classifier"],
                                 tax_DB = config["tax_DB_name"],
                                 rarefaction_value = get_rarefaction_key(config["rarefaction_value"]),
                                 grouping_column = config["grouping_column"][0]))


basic_plots.append(expand("{denoiser}/5_visualization/{classifier}_{tax_DB}/rarefaction_curve/{filter_lineage}_in_{filter_tax_rank}/{rarefaction_value}_rarefaction_curve_{grouping_column}.png",
                                 denoiser = config["denoiser"],
                                 classifier = config["classifier"],
                                 tax_DB = config["tax_DB_name"],
                                 rarefaction_value = get_rarefaction_key(config["rarefaction_value"]),
                                 filter_tax_rank = config["filter_tax_rank"],
                                 filter_lineage = config["filter_lineage"],
                                 grouping_column = config["grouping_column"][0]))



basic_plots.append(expand("{denoiser}/5_visualization/{classifier}_{tax_DB}/KRONA/{grouping_key}.html",
                                 denoiser = config["denoiser"],
                                 classifier = config["classifier"],
                                 tax_DB = config["tax_DB_name"],
                                 grouping_key = get_grouping_key(config["grouping_column"])))



### Advances Phyloseq set of output
#### Not filtered
phyloseq = expand("{denoiser}/4_physeq/{classifier}_{tax_DB}/nofiltering/{rarefaction_value}/{collapse_key}/base_with_tree.rds",
                   denoiser = config["denoiser"],
                   classifier = config["classifier"],
                   tax_DB = config["tax_DB_name"],
                   collapse_key = get_taxa_collapse_level_key(config["phyloseq_tax_ranks"]),
                   rarefaction_value = get_rarefaction_key(config["rarefaction_value"]))


phyloseq.append(expand("{denoiser}/4_physeq/{classifier}_{tax_DB}/nofiltering/{rarefaction_value}/{collapse_key}/base_with_tree.rds",
                   denoiser = config["denoiser"],
                   classifier = config["classifier"],
                   tax_DB = config["tax_DB_name"],
                   collapse_key = get_taxa_collapse_level_key(config["phyloseq_tax_ranks"]),
                   rarefaction_value = get_rarefaction_key(config["rarefaction_value"])))

phyloseq.append(expand("{denoiser}/4_physeq/{classifier}_{tax_DB}/nofiltering/{rarefaction_value}/{collapse_key}/base_with_tree_norm{norm_value}_abund{abund_value}_prev{prev_value}.rds",
                   denoiser = config["denoiser"],
                   classifier = config["classifier"],
                   tax_DB = config["tax_DB_name"],
                   collapse_key = get_taxa_collapse_level_key(config["phyloseq_tax_ranks"]),
                   rarefaction_value = get_rarefaction_key(config["rarefaction_value"]),
                   norm_value = config["normalization"], 
                   abund_value = config["min_abundance"], 
                   prev_value = config["min_prevalence"]))

#### Tax filtered
phyloseq.append(expand("{denoiser}/4_physeq/{classifier}_{tax_DB}/{filter_lineage}_in_{filter_tax_rank}/{rarefaction_value}/{collapse_key}/base_with_tree.rds",
                   denoiser = config["denoiser"],
                   classifier = config["classifier"],
                   tax_DB = config["tax_DB_name"],
                   collapse_key = get_taxa_collapse_level_key(config["phyloseq_tax_ranks"]),
                   rarefaction_value = get_rarefaction_key(config["rarefaction_value"]),
                   filter_tax_rank = config["filter_tax_rank"],
                   filter_lineage = config["filter_lineage"]))


phyloseq.append(expand("{denoiser}/4_physeq/{classifier}_{tax_DB}/{filter_lineage}_in_{filter_tax_rank}/{rarefaction_value}/{collapse_key}/base_with_tree.rds",
                   denoiser = config["denoiser"],
                   classifier = config["classifier"],
                   tax_DB = config["tax_DB_name"],
                   collapse_key = get_taxa_collapse_level_key(config["phyloseq_tax_ranks"]),
                   rarefaction_value = get_rarefaction_key(config["rarefaction_value"]),
                   filter_tax_rank = config["filter_tax_rank"],
                   filter_lineage = config["filter_lineage"]))

phyloseq.append(expand("{denoiser}/4_physeq/{classifier}_{tax_DB}/{filter_lineage}_in_{filter_tax_rank}/{rarefaction_value}/{collapse_key}/base_with_tree_norm{norm_value}_abund{abund_value}_prev{prev_value}.rds",
                   denoiser = config["denoiser"],
                   classifier = config["classifier"],
                   tax_DB = config["tax_DB_name"],
                   collapse_key = get_taxa_collapse_level_key(config["phyloseq_tax_ranks"]),
                   rarefaction_value = get_rarefaction_key(config["rarefaction_value"]),
                   filter_tax_rank = config["filter_tax_rank"],
                   filter_lineage = config["filter_lineage"],
                   norm_value = config["normalization"], 
                   abund_value = config["min_abundance"], 
                   prev_value = config["min_prevalence"]))


### Potentially spurious ASVs
phyloseq.append(expand("{denoiser}/4_physeq/{classifier}_{tax_DB}/{filter_lineage}_in_{filter_tax_rank}/norarefaction/no_collapse/base_export/tree_treeshrink/output.txt",
            denoiser = config["denoiser"],
            classifier = config["classifier"],
            tax_DB = config["tax_DB_name"],
            collapse_key = get_taxa_collapse_level_key(config["phyloseq_tax_ranks"]),
            rarefaction_value = get_rarefaction_key(config["rarefaction_value"]),
            filter_tax_rank = config["filter_tax_rank"],
            filter_lineage = config["filter_lineage"],
            norm_value = config["normalization"], 
            abund_value = config["min_abundance"], 
            prev_value = config["min_prevalence"]))


#### Melted Phyloseq set of output
#### Not filtered
phyloseq_melted = expand("{denoiser}/4_physeq/{classifier}_{tax_DB}/nofiltering/{rarefaction_value}/{collapse_key}/base_with_tree_melted.tsv",
                   denoiser = config["denoiser"],
                   classifier = config["classifier"],
                   tax_DB = config["tax_DB_name"],
                   collapse_key = get_taxa_collapse_level_key(config["phyloseq_tax_ranks"]),
                   rarefaction_value = get_rarefaction_key(config["rarefaction_value"]))


phyloseq_melted.append(expand("{denoiser}/4_physeq/{classifier}_{tax_DB}/nofiltering/{rarefaction_value}/{collapse_key}/base_with_tree_melted.tsv",
                   denoiser = config["denoiser"],
                   classifier = config["classifier"],
                   tax_DB = config["tax_DB_name"],
                   collapse_key = get_taxa_collapse_level_key(config["phyloseq_tax_ranks"]),
                   rarefaction_value = get_rarefaction_key(config["rarefaction_value"])))

phyloseq_melted.append(expand("{denoiser}/4_physeq/{classifier}_{tax_DB}/nofiltering/{rarefaction_value}/{collapse_key}/base_with_tree_norm{norm_value}_abund{abund_value}_prev{prev_value}_melted.tsv",
                   denoiser = config["denoiser"],
                   classifier = config["classifier"],
                   tax_DB = config["tax_DB_name"],
                   collapse_key = get_taxa_collapse_level_key(config["phyloseq_tax_ranks"]),
                   rarefaction_value = get_rarefaction_key(config["rarefaction_value"]),
                   norm_value = config["normalization"], 
                   abund_value = config["min_abundance"], 
                   prev_value = config["min_prevalence"]))

#### Tax filtered
phyloseq_melted.append(expand("{denoiser}/4_physeq/{classifier}_{tax_DB}/{filter_lineage}_in_{filter_tax_rank}/{rarefaction_value}/{collapse_key}/base_with_tree_melted.tsv",
                   denoiser = config["denoiser"],
                   classifier = config["classifier"],
                   tax_DB = config["tax_DB_name"],
                   collapse_key = get_taxa_collapse_level_key(config["phyloseq_tax_ranks"]),
                   rarefaction_value = get_rarefaction_key(config["rarefaction_value"]),
                   filter_tax_rank = config["filter_tax_rank"],
                   filter_lineage = config["filter_lineage"]))


phyloseq_melted.append(expand("{denoiser}/4_physeq/{classifier}_{tax_DB}/{filter_lineage}_in_{filter_tax_rank}/{rarefaction_value}/{collapse_key}/base_with_tree_melted.tsv",
                   denoiser = config["denoiser"],
                   classifier = config["classifier"],
                   tax_DB = config["tax_DB_name"],
                   collapse_key = get_taxa_collapse_level_key(config["phyloseq_tax_ranks"]),
                   rarefaction_value = get_rarefaction_key(config["rarefaction_value"]),
                   filter_tax_rank = config["filter_tax_rank"],
                   filter_lineage = config["filter_lineage"]))

phyloseq_melted.append(expand("{denoiser}/4_physeq/{classifier}_{tax_DB}/{filter_lineage}_in_{filter_tax_rank}/{rarefaction_value}/{collapse_key}/base_with_tree_norm{norm_value}_abund{abund_value}_prev{prev_value}_melted.tsv",
                   denoiser = config["denoiser"],
                   classifier = config["classifier"],
                   tax_DB = config["tax_DB_name"],
                   collapse_key = get_taxa_collapse_level_key(config["phyloseq_tax_ranks"]),
                   rarefaction_value = get_rarefaction_key(config["rarefaction_value"]),
                   filter_tax_rank = config["filter_tax_rank"],
                   filter_lineage = config["filter_lineage"],
                   norm_value = config["normalization"], 
                   abund_value = config["min_abundance"], 
                   prev_value = config["min_prevalence"]))


#### Transposed output tables
#### Count-table transposed format - filtered
transposed_output = expand("{denoiser}/4_physeq/{classifier}_{tax_DB}/{filter_lineage}_in_{filter_tax_rank}/{rarefaction_value}/{collapse_key}/base_with_tree_norm{norm_value}_abund{abund_value}_prev{prev_value}_export/count_table_transposed.txt",
            denoiser = config["denoiser"],
            classifier = config["classifier"],
            tax_DB = config["tax_DB_name"],
            collapse_key = get_taxa_collapse_level_key(config["phyloseq_tax_ranks"]),
            rarefaction_value = get_rarefaction_key(config["rarefaction_value"]),
            filter_tax_rank = config["filter_tax_rank"],
            filter_lineage = config["filter_lineage"],
            norm_value = config["normalization"], 
            abund_value = config["min_abundance"], 
            prev_value = config["min_prevalence"])


#### Count-table transposed format - not filtered
transposed_output.append(expand("{denoiser}/4_physeq/{classifier}_{tax_DB}/nofiltering/{rarefaction_value}/{collapse_key}/base_with_tree_norm{norm_value}_abund{abund_value}_prev{prev_value}_export/count_table_transposed.txt",
            denoiser = config["denoiser"],
            classifier = config["classifier"],
            tax_DB = config["tax_DB_name"],
            collapse_key = get_taxa_collapse_level_key(config["phyloseq_tax_ranks"]),
            rarefaction_value = get_rarefaction_key(config["rarefaction_value"]),
            filter_tax_rank = config["filter_tax_rank"],
            filter_lineage = config["filter_lineage"],
            norm_value = config["normalization"], 
            abund_value = config["min_abundance"], 
            prev_value = config["min_prevalence"]))


### Qiime2 output
#### Qiime2 interactive visualization
Qiime2_vis_qzv = expand("{denoiser}/2_denoised/dna-sequences.fasta",
                        denoiser = config["denoiser"])
Qiime2_vis_qzv.append(expand("{denoiser}/2_denoised/count-table.qzv",
                             denoiser = config["denoiser"]))
Qiime2_vis_qzv.append(expand("{denoiser}/2_denoised/rep-seqs.qzv",
                             denoiser = config["denoiser"]))
Qiime2_vis_qzv.append(expand("{denoiser}/3_classified/{classifier}_{tax_DB}/dna-sequences_tax_assignments.qzv",
                      denoiser = config["denoiser"],
                      classifier = config["classifier"],
                      tax_DB = config["tax_DB_name"]))


### Picrust2
picrust2 = expand("{denoiser}/6_picrust2/{classifier}_{tax_DB}/{filter_lineage}_in_{filter_tax_rank}/{rarefaction_value}/picrust/",
                    denoiser = config["denoiser"],
                    classifier = config["classifier"],
                    tax_DB = config["tax_DB_name"],
                    filter_tax_rank = config["filter_tax_rank"],
                    filter_lineage = config["filter_lineage"],
                    rarefaction_value = get_rarefaction_key(config["rarefaction_value"])
                    )

