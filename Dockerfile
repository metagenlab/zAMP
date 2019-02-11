FROM continuumio/miniconda3:4.5.12

RUN conda config --add channels defaults && conda config --add channels conda-forge && conda config --add channels bioconda

# RUN useradd -r -u 1080 pipeline_user

# Install most of needed packages
RUN conda create -n r_visualization

RUN source activate r_visualization

RUN conda env update -f=envs/r_visualization.yml

# randomcoloR dependancy
RUN export DEBIAN_FRONTEND=noninteractive TERM=linux && \
  apt-get update && \
  apt-get -y --no-install-recommends install libv8-dev

RUN source activate r_visualization

# randomcoloR
RUN echo 'install.packages(randomcoloR, repos="http://cran.us.r-project.org", dependencies=TRUE)' > /tmp/packages.R \ && Rscript /tmp/packages.R
