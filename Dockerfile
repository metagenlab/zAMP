FROM continuumio/miniconda3:4.5.12

ADD envs/r_visualization.yml /tmp/r_visualization.yml
RUN conda env create -f /tmp/r_visualization.yml -n r_visualization

# Pull the environment name out of the r_visualization.yml
RUN echo "source activate r_visualization" > ~/.bashrc
ENV PATH /opt/conda/envs/r_visualization/bin:$PATH

RUN source activate r_visualization

RUN export DEBIAN_FRONTEND=noninteractive TERM=linux && \
  apt-get update && \
  apt-get update && \
  apt-get -y --no-install-recommends install libv8-dev

RUN echo 'install.packages(randomcoloR, repos="http://cran.us.r-project.org", dependencies=TRUE)' > /tmp/packages.R \ && Rscript /tmp/packages.R
