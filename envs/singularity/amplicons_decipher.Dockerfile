FROM continuumio/miniconda3:4.7.12

################## METADATA ######################
 
LABEL base.image="miniconda3:4.7.12"
LABEL version="v.1.0"
LABEL software="DECIPHER"
LABEL software.version="2.12"
LABEL description="DECIPHER (https://doi.org/10.1186/s40168-018-0521-5) with dplyr for output formatting"
LABEL tags="Genomics"
 
################## MAINTAINER ######################
 
MAINTAINER Valentin Scherz
 
################## INSTALLATION ######################
ENV DEBIAN_FRONTEND noninteractive

COPY ./envs/decipher.yml ./decipher.yml

RUN conda config --add channels bioconda && \
    conda config --add channels conda-forge && \ 
    conda update conda && \
    conda env create -f decipher.yml && \
    conda clean --all --yes


RUN conda init bash
ENTRYPOINT ["/bin/bash"]
ENV PATH /opt/conda/envs/decipher/bin:$PATH