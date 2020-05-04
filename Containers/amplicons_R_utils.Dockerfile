FROM continuumio/miniconda:4.7.12


########################### Install with conda some basic usefull R packages ##############################

COPY ./envs/amplicons_R_utils.yml /tmp/conda_definition.yml



RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda config --add channels biocore && \
    conda env update --name base --file /tmp/conda_definition.yml

RUN conda init bash

ENTRYPOINT ["/bin/bash"]
