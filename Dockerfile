FROM continuumio/miniconda3:4.5.12

ADD envs/r_visualization.yml /tmp/r_visualization.yml
RUN conda env create -f /tmp/r_visualization.yml

# Pull the environment name out of the r_visualization.yml
RUN echo "source activate r_visualization" > ~/.bashrc
ENV PATH /opt/conda/envs/r_visualization/bin:$PATH
