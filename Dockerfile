FROM continuumio/miniconda3:4.5.12

RUN conda create -n r_visualization python=3.6
RUN echo "source activate r_visualization" > ~/.bashrc
ENV PATH /opt/conda/envs/env/bin:$PATH

RUN conda config --add channels defaults && conda config --add channels conda-forge && conda config --add channels bioconda

# randomcoloR dependancy
RUN export DEBIAN_FRONTEND=noninteractive TERM=linux && \
  apt-get update && \
  apt-get -y --no-install-recommends install libv8-dev

#RUN source activate r_visualization

# randomcoloR
RUN echo 'install.packages(randomcoloR, repos="http://cran.us.r-project.org", dependencies=TRUE)' > /tmp/packages.R \ && Rscript /tmp/packages.R


# RUN useradd -r -u 1080 pipeline_user
