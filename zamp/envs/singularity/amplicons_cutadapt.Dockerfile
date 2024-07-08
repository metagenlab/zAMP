FROM continuumio/miniconda3:4.7.12


################## METADATA ######################
 
LABEL base.image="miniconda3:4.7.12"
LABEL version="v.3.0"
LABEL software="amplicons_cutadapt"
LABEL software.version="4.7"
LABEL description="Cutadapt http://dx.doi.org/10.14806/ej.17.1.200 with Biopython"
LABEL tags="Genomics"
 
################## MAINTAINER ######################
 
MAINTAINER Valentin Scherz
 
################## INSTALLATION ######################
ENV DEBIAN_FRONTEND noninteractive


COPY ./envs/cutadapt.yml ./cutadapt.yml

RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \ 
    conda update conda && \
    conda env create -f cutadapt.yml && \
    conda clean --all --yes


RUN conda init bash
ENTRYPOINT ["/bin/bash"]
ENV PATH /opt/conda/envs/cutadapt/bin:$PATH
