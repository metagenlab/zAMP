FROM continuumio/miniconda:4.6.14

## Install Qiime1 and rdp
RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda config --add channels hcc && \
    conda install blast=2.9.0 java-jdk conda=4.6.14 qiime=1.9.1=np110py27_1 rdp-classifier=2.2



