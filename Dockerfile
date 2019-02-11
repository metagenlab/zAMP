FROM continuumio/miniconda3:4.5.12

ADD envs/r_visualization.yml /tmp/r_visualization.yml
RUN conda env create -f /tmp/r_visualization.yml

# Pull the environment name out of the r_visualization.yml
RUN echo "source activate $(head -1 /tmp/r_visualization.yml | cut -d' ' -f2)" > ~/.bashrc
ENV PATH /opt/conda/envs/$(head -1 /tmp/r_visualization.yml | cut -d' ' -f2)/bin:$PATH
