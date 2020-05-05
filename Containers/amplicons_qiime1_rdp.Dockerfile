FROM continuumio/miniconda:4.6.14

################## METADATA ######################
 
LABEL base.image="miniconda:4.6.14"
LABEL version="v.1.0"
LABEL software="QIIME_RDP"
LABEL software.version="1.9.1/2.2"
LABEL description="Qiime (https://dx.doi.org/10.1038%2Fnmeth.f.303) and RDP (https://dx.doi.org/10.1128%2FAEM.00062-07)"
LABEL tags="Genomics"
 
################## MAINTAINER ######################
 
MAINTAINER Valentin Scherz
 
################## INSTALLATION ######################
ENV DEBIAN_FRONTEND noninteractive


COPY ./envs/QIIME1.yml ./QIIME1.yml

RUN conda config --add channels conda-forge && \
    conda config --add channels bioconda && \ 
    conda config --add channels hcc && \ 
    conda env create -f QIIME1.yml && \
    conda clean --all --yes


RUN conda init bash
ENTRYPOINT ["/bin/bash"]
ENV PATH /opt/conda/envs/QIIME1/bin:$PATH
