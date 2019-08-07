
ARG VERSION
ARG TEST_CPU

FROM metagenlab/amplicon_pipeline:$VERSION

#################### Run the pipeline to test it #####################
## Recover input files
RUN cp ${pipeline_folder}/data/validation_datasets/config.yml ./ && \
    cp ${pipeline_folder}/data/validation_datasets/config_in_silico.yml ./ && \
    cp ${pipeline_folder}/data/validation_datasets/input_table_insilico.tsv ./ && \
    cp ${pipeline_folder}/data/validation_datasets/validation_set.tsv ./

## Run the pipeline to test it
RUN snakemake --snakefile ${pipeline_folder}/Snakefile --cores $TEST_CPU --resources max_copy=4 --use-conda --conda-prefix /opt/conda/ --configfile config.yml all PICRUSt2_output
