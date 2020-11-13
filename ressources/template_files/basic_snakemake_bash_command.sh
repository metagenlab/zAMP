snakemake \
    --snakefile <path/to/pipeline>/microbiome16S_pipeline>/Snakefile \
    --use-singularity --singularity-prefix <path/to/singularity/containers/storage/> \
    #or --use-conda --conda-prefix <path/to/conda/environements/storage/> \
    --cores <number_of_CPU_allocated_to_the_analysis> \
    --configfile config.yaml \
    --resources max_copy=<number_of_file_copied_simultaneously_4_is_a_good_default_value> mem_mb=<available_RAM_memory_in_MB>\
    all #create_filtered_reads_multiqc_report PICRUSt2_output