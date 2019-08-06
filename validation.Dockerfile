
ARG VERSION
ARG TEST_CPU

FROM metagenlab/amplicon_pipeline:$VERSION

#################### Run the pipeline to test it #####################
## Recovert input files
RUN cp ${pipeline_folder}/data/validation_datasets/config.yml ./ && \
    cp ${pipeline_folder}/data/validation_datasets/config_in_silico.yml ./ && \
    cp ${pipeline_folder}/data/validation_datasets/input_table_insilico.tsv ./ && \
    cp ${pipeline_folder}/data/validation_datasets/validation_set.tsv ./

## Test the insilico validation
RUN snakemake --snakefile ${pipeline_folder}/Snakefile_validation --cores $TEST_CPU --resources ncbi_requests=2 --use-conda --conda-prefix /opt/conda/ --configfile config_in_silico.yml insilico_validation

## Run the pipeline to test it
RUN snakemake --snakefile ${pipeline_folder}/Snakefile --cores $TEST_CPU --resources max_copy=4 --use-conda --conda-prefix /opt/conda/ --configfile config.yml all PICRUSt2_output
