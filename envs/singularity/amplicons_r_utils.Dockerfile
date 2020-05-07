FROM continuumio/miniconda3:4.7.12

################## METADATA ######################
 
LABEL base.image="miniconda3:4.7.12"
LABEL version="v.1.0"
LABEL software="amplicons_R_utils"
LABEL software.version="1.0"
LABEL description="Various R packages for our amplicon-based metagenomics pipeline"
LABEL tags="Genomics"
 
################## MAINTAINER ######################
 
MAINTAINER Valentin Scherz
 
################## INSTALLATION ######################
ENV DEBIAN_FRONTEND noninteractive

COPY ./envs/amplicons_r_utils.yml ./amplicons_r_utils.yml

RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \ 
    conda update conda && \
    conda env create -f amplicons_r_utils.yml && \
    conda clean --all --yes


RUN conda init bash
ENTRYPOINT ["/bin/bash"]
ENV PATH /opt/conda/envs/amplicons_r_utils/bin:$PATH


