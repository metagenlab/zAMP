snakemake \
    --snakefile <path/to/pipeline>/microbiome16S_pipeline>/DBprocess.Snakefile \
    --use-singularity --singularity-prefix <path/to/singularity/containers/storage/> \
    #or --use-conda --conda-prefix <path/to/conda/environements/storage/> \
    --cores <number_of_CPU_allocated_to_the_analysis> \
    --configfile config_DB.yaml