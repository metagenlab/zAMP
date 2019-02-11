FROM continuumio/miniconda3:4.5.12

ADD envs/r_visualization.yml /tmp/r_visualization.yml
RUN conda env create -f /tmp/r_visualization.yml -n r_visualization

RUN /bin/bash -c 'source activate /opt/conda/envs/r_visualization/'

RUN apt-get update && apt-get install -y --no-install-recommends apt-utils

RUN apt-get update && apt-get -y install libv8-dev && apt-get -y install g++

RUN dpkg -L libv8-dev


RUN /bin/bash -c 'source activate /opt/conda/envs/r_visualization/ && R CMD install --configure-vars='INCLUDE_DIR="/usr/include"',  randomcoloR
