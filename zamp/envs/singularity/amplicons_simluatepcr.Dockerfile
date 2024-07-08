FROM continuumio/miniconda3:4.7.12

################## METADATA ######################
 
LABEL base.image="miniconda3:4.7.12"
LABEL version="v.1.0"
LABEL software="simulatePCR"
LABEL software.version="v.0.9.9"
LABEL description="Patched version of simulate_PCR (https://doi.org/10.1186/1471-2105-15-237)"
LABEL tags="Genomics"
 
################## MAINTAINER ######################
 
MAINTAINER Valentin Scherz
 
################## INSTALLATION ######################

## Install simulatePCR dependancies via Conda
RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda install blast=2.9.0 perl-lwp-simple perl-bioperl java-jdk conda=4.6.14 && \
    conda clean --all --yes

## Set in path a patched version of simulate PCR, DOI: 10.1186/1471-2105-15-237 for amplicons validation
RUN wget --quiet https://github.com/metagenlab/updated_simulate_PCR/archive/v0.9.9.tar.gz -O simulate_PCR.tar.gz && \
    mkdir /opt/simulate_PCR && \
    tar xzf simulate_PCR.tar.gz -C /opt/simulate_PCR &&  \
    mv /opt/simulate_PCR/updated_simulate_PCR-0.9.9/code/simulate_PCR /opt/simulate_PCR && \
    rm simulate_PCR.tar.gz && \
    rm -R /opt/simulate_PCR/updated_simulate_PCR-0.9.9 && \
    chmod +x /opt/simulate_PCR/simulate_PCR

ENV PATH="/opt/simulate_PCR:${PATH}"
ENV PERL5LIB="/opt/conda/lib/site_perl/5.26.2"
ENV PATH="/opt/simulate_PCR:${PATH}"

ENTRYPOINT ["/bin/bash"]
