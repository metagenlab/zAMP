FROM continuumio/miniconda3:4.7.12

################## METADATA ######################
 
LABEL base.image="miniconda3:4.7.12"
LABEL version="v.1.0"
LABEL software="TreeShrink"
LABEL software.version="1.3.7"
LABEL description="TreeShrink to identify outliers in phylogenetic tree. https://github.com/uym2/TreeShrink & https://doi.org/10.1186/s12864-018-4620-2"
LABEL tags="Genomics"
 
################## MAINTAINER ######################
 
MAINTAINER Valentin Scherz
 
################## INSTALLATION ######################
ENV DEBIAN_FRONTEND noninteractive

COPY ./envs/TreeShrink.yml ./TreeShrink.yml

RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \ 
    conda update conda && \
    conda env create -f TreeShrink.yml && \
    conda clean --all --yes


RUN conda init bash
ENTRYPOINT ["/bin/bash"]
ENV PATH /opt/conda/envs/treeshrink/bin:$PATH


