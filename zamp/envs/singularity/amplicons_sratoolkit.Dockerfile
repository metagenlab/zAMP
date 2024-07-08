FROM continuumio/miniconda3:4.7.12


################## METADATA ######################
 
LABEL base.image="miniconda3:4.7.12"
LABEL version="v.1.0"
LABEL software="SRA toolkit"
LABEL software.version="3.10.0"
LABEL description=" SRA toolkit"
LABEL tags="Genomics"
 
################## MAINTAINER ######################
 
MAINTAINER Valentin Scherz
 
################## INSTALLATION ######################
ENV DEBIAN_FRONTEND noninteractive


COPY ./envs/sra-tools.yml ./sra-tools.yml

RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \ 
    conda update conda && \
    conda env create -f sra-tools.yml && \
    conda clean --all --yes


RUN conda init bash
ENTRYPOINT ["/bin/bash"]
ENV PATH /opt/conda/envs/sra-tools/bin:$PATH