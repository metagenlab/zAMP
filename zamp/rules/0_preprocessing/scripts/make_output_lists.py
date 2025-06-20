## Define lits of outputs as that will be called by function from the Snakefile depending from the values in config.

### Import list of functions to handle the output
include: "make_output_fcts.py"

### Basic output
#### MuliQC
MultiQC = expand(
    os.path.join("{qc}", "{run}_multiqc_raw_reads_report.html"),
    qc=dir.out.qc,
    run=set(SAMPLES["run"]),
)
MultiQC.append(f"{dir.out.qc}/multiqc_raw_reads_report.html")


#### Light pipeline output
light_output = expand(
    os.path.join(dir.out.base, "{denoiser}", "2_denoised", "dna-sequences.fasta"),
    denoiser=DENOISER,
)
light_output.append(
    expand(
        os.path.join(dir.out.base, "{denoiser}", "2_denoised", "count_table.tsv"),
        denoiser=DENOISER,
    )
)
light_output.append(
    expand(
        os.path.join(
            dir.out.base,
            "{denoiser}",
            "3_classified",
            "{classifier}_{tax_DB}",
            "taxonomy.qzv",
        ),
        classifier=CLASSIFIER,
        denoiser=DENOISER,
        tax_DB=DBNAME,
    )
)

#### DADA2 stat table
DADA2_stats_tables = os.path.join(
    dir.out.base, "DADA2", "2_denoised", "DADA2_denoising_stats.tsv"
)


### Basic evaluation plots
basic_plots = expand(
    os.path.join(
        dir.out.base,
        "{denoiser}",
        "5_visualization",
        "{classifier}_{tax_DB}",
        "reads",
        "reads_plot_with_filtered.png",
    ),
    classifier=CLASSIFIER,
    denoiser=DENOISER,
    tax_DB=DBNAME,
)


basic_plots.append(
    expand(
        os.path.join(
            dir.out.base,
            "{denoiser}",
            "5_visualization",
            "{classifier}_{tax_DB}",
            "reads",
            "reads_plot_with_{filter_lineage}_in_{filter_tax_rank}_filtered.png",
        ),
        denoiser=DENOISER,
        classifier=CLASSIFIER,
        tax_DB=DBNAME,
        filter_tax_rank=KEEP_RANK,
        filter_lineage=KEEP_TAXA,
    )
)


basic_plots.append(
    expand(
        os.path.join(
            dir.out.base,
            "{denoiser}",
            "5_visualization",
            "{classifier}_{tax_DB}",
            "rarefaction_curve",
            "nofiltering",
            "{rarefaction_value}_rarefaction_curve_{grouping_column}.png",
        ),
        denoiser=DENOISER,
        classifier=CLASSIFIER,
        tax_DB=DBNAME,
        rarefaction_value=get_rarefaction_key(RAREFACTION),
        grouping_column="sample_group",
    )
)


basic_plots.append(
    expand(
        os.path.join(
            dir.out.base,
            "{denoiser}",
            "5_visualization",
            "{classifier}_{tax_DB}",
            "rarefaction_curve",
            "{filter_lineage}_in_{filter_tax_rank}",
            "{rarefaction_value}_rarefaction_curve_{grouping_column}.png",
        ),
        denoiser=DENOISER,
        classifier=CLASSIFIER,
        tax_DB=DBNAME,
        rarefaction_value=get_rarefaction_key(RAREFACTION),
        filter_tax_rank=KEEP_RANK,
        filter_lineage=KEEP_TAXA,
        grouping_column="sample_group",
    )
)


basic_plots.append(
    expand(
        os.path.join(
            dir.out.base,
            "{denoiser}",
            "5_visualization",
            "{classifier}_{tax_DB}",
            "KRONA",
            "{grouping_key}.html",
        ),
        denoiser=DENOISER,
        classifier=CLASSIFIER,
        tax_DB=DBNAME,
        grouping_key=get_grouping_key("sample_group"),
    )
)


### Advances Phyloseq set of output
#### Not filtered

phyloseq = expand(
    os.path.join(
        dir.out.base,
        "{denoiser}",
        "4_physeq",
        "{classifier}_{tax_DB}",
        "nofiltering",
        "{rarefaction_value}",
        "{collapse_key}",
        "base_with_tree.rds",
    ),
    denoiser=DENOISER,
    classifier=CLASSIFIER,
    tax_DB=DBNAME,
    collapse_key=get_taxa_collapse_level_key(PHYSEQ_RANK),
    rarefaction_value=get_rarefaction_key(RAREFACTION),
)


phyloseq.append(
    expand(
        os.path.join(
            dir.out.base,
            "{denoiser}",
            "4_physeq",
            "{classifier}_{tax_DB}",
            "nofiltering",
            "{rarefaction_value}",
            "{collapse_key}",
            "base_with_tree.rds",
        ),
        denoiser=DENOISER,
        classifier=CLASSIFIER,
        tax_DB=DBNAME,
        collapse_key=get_taxa_collapse_level_key(PHYSEQ_RANK),
        rarefaction_value=get_rarefaction_key(RAREFACTION),
    )
)

phyloseq.append(
    expand(
        os.path.join(
            dir.out.base,
            "{denoiser}",
            "4_physeq",
            "{classifier}_{tax_DB}",
            "nofiltering",
            "{rarefaction_value}",
            "{collapse_key}",
            "base_with_tree_norm{norm_value}_abund{abund_value}_prev{prev_value}.rds",
        ),
        denoiser=DENOISER,
        classifier=CLASSIFIER,
        tax_DB=DBNAME,
        collapse_key=get_taxa_collapse_level_key(PHYSEQ_RANK),
        rarefaction_value=get_rarefaction_key(RAREFACTION),
        norm_value=NORM,
        abund_value=MIN_COUNT,
        prev_value=MIN_PREV,
    )
)

#### Tax filtered
phyloseq.append(
    expand(
        os.path.join(
            dir.out.base,
            "{denoiser}",
            "4_physeq",
            "{classifier}_{tax_DB}",
            "{filter_lineage}_in_{filter_tax_rank}",
            "{rarefaction_value}",
            "{collapse_key}",
            "base_with_tree.rds",
        ),
        denoiser=DENOISER,
        classifier=CLASSIFIER,
        tax_DB=DBNAME,
        collapse_key=get_taxa_collapse_level_key(PHYSEQ_RANK),
        rarefaction_value=get_rarefaction_key(RAREFACTION),
        filter_tax_rank=KEEP_RANK,
        filter_lineage=KEEP_TAXA,
    )
)


phyloseq.append(
    expand(
        os.path.join(
            dir.out.base,
            "{denoiser}",
            "4_physeq",
            "{classifier}_{tax_DB}",
            "{filter_lineage}_in_{filter_tax_rank}",
            "{rarefaction_value}",
            "{collapse_key}",
            "base_with_tree.rds",
        ),
        denoiser=DENOISER,
        classifier=CLASSIFIER,
        tax_DB=DBNAME,
        collapse_key=get_taxa_collapse_level_key(PHYSEQ_RANK),
        rarefaction_value=get_rarefaction_key(RAREFACTION),
        filter_tax_rank=KEEP_RANK,
        filter_lineage=KEEP_TAXA,
    )
)

phyloseq.append(
    expand(
        os.path.join(
            dir.out.base,
            "{denoiser}",
            "4_physeq",
            "{classifier}_{tax_DB}",
            "{filter_lineage}_in_{filter_tax_rank}",
            "{rarefaction_value}",
            "{collapse_key}",
            "base_with_tree_norm{norm_value}_abund{abund_value}_prev{prev_value}.rds",
        ),
        denoiser=DENOISER,
        classifier=CLASSIFIER,
        tax_DB=DBNAME,
        collapse_key=get_taxa_collapse_level_key(PHYSEQ_RANK),
        rarefaction_value=get_rarefaction_key(RAREFACTION),
        filter_tax_rank=KEEP_RANK,
        filter_lineage=KEEP_TAXA,
        norm_value=NORM,
        abund_value=MIN_COUNT,
        prev_value=MIN_PREV,
    )
)


### Potentially spurious ASVs
phyloseq.append(
    expand(
        os.path.join(
            dir.out.base,
            "{denoiser}",
            "4_physeq",
            "{classifier}_{tax_DB}",
            "{filter_lineage}_in_{filter_tax_rank}",
            "norarefaction",
            "no_collapse",
            "base_export",
            "tree_treeshrink",
            "output.txt",
        ),
        denoiser=DENOISER,
        classifier=CLASSIFIER,
        tax_DB=DBNAME,
        collapse_key=get_taxa_collapse_level_key(PHYSEQ_RANK),
        rarefaction_value=get_rarefaction_key(RAREFACTION),
        filter_tax_rank=KEEP_RANK,
        filter_lineage=KEEP_TAXA,
        norm_value=NORM,
        abund_value=MIN_COUNT,
        prev_value=MIN_PREV,
    )
)


#### Melted Phyloseq set of output
#### Not filtered
phyloseq_melted = expand(
    os.path.join(
        dir.out.base,
        "{denoiser}",
        "4_physeq",
        "{classifier}_{tax_DB}",
        "nofiltering",
        "{rarefaction_value}",
        "{collapse_key}",
        "base_with_tree_melted.tsv",
    ),
    denoiser=DENOISER,
    classifier=CLASSIFIER,
    tax_DB=DBNAME,
    collapse_key=get_taxa_collapse_level_key(PHYSEQ_RANK),
    rarefaction_value=get_rarefaction_key(RAREFACTION),
)


phyloseq_melted.append(
    expand(
        os.path.join(
            dir.out.base,
            "{denoiser}",
            "4_physeq",
            "{classifier}_{tax_DB}",
            "nofiltering",
            "{rarefaction_value}",
            "{collapse_key}",
            "base_with_tree_melted.tsv",
        ),
        denoiser=DENOISER,
        classifier=CLASSIFIER,
        tax_DB=DBNAME,
        collapse_key=get_taxa_collapse_level_key(PHYSEQ_RANK),
        rarefaction_value=get_rarefaction_key(RAREFACTION),
    )
)

phyloseq_melted.append(
    expand(
        os.path.join(
            dir.out.base,
            "{denoiser}",
            "4_physeq",
            "{classifier}_{tax_DB}",
            "nofiltering",
            "{rarefaction_value}",
            "{collapse_key}",
            "base_with_tree_norm{norm_value}_abund{abund_value}_prev{prev_value}_melted.tsv",
        ),
        denoiser=DENOISER,
        classifier=CLASSIFIER,
        tax_DB=DBNAME,
        collapse_key=get_taxa_collapse_level_key(PHYSEQ_RANK),
        rarefaction_value=get_rarefaction_key(RAREFACTION),
        norm_value=NORM,
        abund_value=MIN_COUNT,
        prev_value=MIN_PREV,
    )
)

#### Tax filtered
phyloseq_melted.append(
    expand(
        os.path.join(
            dir.out.base,
            "{denoiser}",
            "4_physeq",
            "{classifier}_{tax_DB}",
            "{filter_lineage}_in_{filter_tax_rank}",
            "{rarefaction_value}",
            "{collapse_key}",
            "base_with_tree_melted.tsv",
        ),
        denoiser=DENOISER,
        classifier=CLASSIFIER,
        tax_DB=DBNAME,
        collapse_key=get_taxa_collapse_level_key(PHYSEQ_RANK),
        rarefaction_value=get_rarefaction_key(RAREFACTION),
        filter_tax_rank=KEEP_RANK,
        filter_lineage=KEEP_TAXA,
    )
)


phyloseq_melted.append(
    expand(
        os.path.join(
            dir.out.base,
            "{denoiser}",
            "4_physeq",
            "{classifier}_{tax_DB}",
            "{filter_lineage}_in_{filter_tax_rank}",
            "{rarefaction_value}",
            "{collapse_key}",
            "base_with_tree_melted.tsv",
        ),
        denoiser=DENOISER,
        classifier=CLASSIFIER,
        tax_DB=DBNAME,
        collapse_key=get_taxa_collapse_level_key(PHYSEQ_RANK),
        rarefaction_value=get_rarefaction_key(RAREFACTION),
        filter_tax_rank=KEEP_RANK,
        filter_lineage=KEEP_TAXA,
    )
)

phyloseq_melted.append(
    expand(
        os.path.join(
            dir.out.base,
            "{denoiser}",
            "4_physeq",
            "{classifier}_{tax_DB}",
            "{filter_lineage}_in_{filter_tax_rank}",
            "{rarefaction_value}",
            "{collapse_key}",
            "base_with_tree_norm{norm_value}_abund{abund_value}_prev{prev_value}_melted.tsv",
        ),
        denoiser=DENOISER,
        classifier=CLASSIFIER,
        tax_DB=DBNAME,
        collapse_key=get_taxa_collapse_level_key(PHYSEQ_RANK),
        rarefaction_value=get_rarefaction_key(RAREFACTION),
        filter_tax_rank=KEEP_RANK,
        filter_lineage=KEEP_TAXA,
        norm_value=NORM,
        abund_value=MIN_COUNT,
        prev_value=MIN_PREV,
    )
)


#### Transposed output tables
#### Count-table transposed format - filtered
transposed_output = expand(
    os.path.join(
        dir.out.base,
        "{denoiser}",
        "4_physeq",
        "{classifier}_{tax_DB}",
        "{filter_lineage}_in_{filter_tax_rank}",
        "{rarefaction_value}",
        "{collapse_key}",
        "base_with_tree_norm{norm_value}_abund{abund_value}_prev{prev_value}_export",
        "count_table_transposed.txt",
    ),
    denoiser=DENOISER,
    classifier=CLASSIFIER,
    tax_DB=DBNAME,
    collapse_key=get_taxa_collapse_level_key(PHYSEQ_RANK),
    rarefaction_value=get_rarefaction_key(RAREFACTION),
    filter_tax_rank=KEEP_RANK,
    filter_lineage=KEEP_TAXA,
    norm_value=NORM,
    abund_value=MIN_COUNT,
    prev_value=MIN_PREV,
)


#### Count-table transposed format - not filtered
transposed_output.append(
    expand(
        os.path.join(
            dir.out.base,
            "{denoiser}",
            "4_physeq",
            "{classifier}_{tax_DB}",
            "{filter_lineage}_in_{filter_tax_rank}",
            "{rarefaction_value}",
            "{collapse_key}",
            "base_with_tree_norm{norm_value}_abund{abund_value}_prev{prev_value}_export",
            "count_table_transposed.txt",
        ),
        denoiser=DENOISER,
        classifier=CLASSIFIER,
        tax_DB=DBNAME,
        collapse_key=get_taxa_collapse_level_key(PHYSEQ_RANK),
        rarefaction_value=get_rarefaction_key(RAREFACTION),
        filter_tax_rank=KEEP_RANK,
        filter_lineage=KEEP_TAXA,
        norm_value=NORM,
        abund_value=MIN_COUNT,
        prev_value=MIN_PREV,
    )
)


### Qiime2 output
#### Qiime2 interactive visualization
Qiime2_vis_qzv = expand(
    os.path.join(dir.out.base, "{denoiser}", "2_denoised", "dna-sequences.fasta"),
    denoiser=DENOISER,
)
Qiime2_vis_qzv.append(
    expand(
        os.path.join(dir.out.base, "{denoiser}", "2_denoised", "count-table.qzv"),
        denoiser=DENOISER,
    )
)
Qiime2_vis_qzv.append(
    expand(
        os.path.join(dir.out.base, "{denoiser}", "2_denoised", "rep-seqs.qzv"),
        denoiser=DENOISER,
    )
)
Qiime2_vis_qzv.append(
    expand(
        os.path.join(
            dir.out.base,
            "{denoiser}",
            "3_classified",
            "{classifier}_{tax_DB}",
            "taxonomy.qzv",
        ),
        denoiser=DENOISER,
        classifier=CLASSIFIER,
        tax_DB=DBNAME,
    )
)


### Picrust2
picrust2 = expand(
    os.path.join(
        dir.out.base,
        "{denoiser}",
        "6_picrust2",
        "{classifier}_{tax_DB}",
        "{filter_lineage}_in_{filter_tax_rank}",
        "{rarefaction_value}",
        "picrust",
    ),
    denoiser=DENOISER,
    classifier=CLASSIFIER,
    tax_DB=DBNAME,
    filter_tax_rank=KEEP_RANK,
    filter_lineage=KEEP_TAXA,
    rarefaction_value=get_rarefaction_key(RAREFACTION),
)
