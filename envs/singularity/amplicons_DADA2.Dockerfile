FROM continuumio/miniconda3:4.7.12

################## METADATA ######################
 
LABEL base.image="miniconda3:4.7.12"
LABEL version="v.2.0"
LABEL software="DADA2"
LABEL software.version="1.26.0"
LABEL description="DADA2 (https://doi.org/10.1038/nmeth.3869) with msfont to fix missing characters in plots"
LABEL tags="Genomics"
 
################## MAINTAINER ######################
 
MAINTAINER Valentin Scherz
 
################## INSTALLATION ######################
ENV DEBIAN_FRONTEND noninteractive

COPY ./envs/DADA2_in_R.yml ./DADA2_in_R.yml

RUN conda config --add channels bioconda && \
    conda config --add channels conda-forge && \ 
    conda update conda && \
    conda env create -f DADA2_in_R.yml && \
    conda clean --all --yes


RUN conda init bash
ENTRYPOINT ["/bin/bash"]
ENV PATH /opt/conda/envs/DADA2/bin:$PATH