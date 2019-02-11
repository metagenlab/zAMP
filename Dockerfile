FROM continuumio/miniconda3:4.5.12

RUN conda config --add channels defaults && conda config --add channels conda-forge && conda config --add channels bioconda

RUN conda info

COPY envs/R_visualization.yml .

RUN cat R_visualization.yml

# RUN useradd -r -u 1080 pipeline_user

# Install most of needed packages
RUN conda create -n R_visualization
RUN echo "source activate R_visualization" > ~/.bashrc
RUN conda install -f ./R_visualization.yml

# Install randomcoloR dependancy
RUN export DEBIAN_FRONTEND=noninteractive TERM=linux && \
  apt-get update && \
  apt-get -y --no-install-recommends install libv8-dev

# Install randomcoloR
RUN echo 'install.packages(randomcoloR, repos="http://cran.us.r-project.org", dependencies=TRUE)' > /tmp/packages.R \ && Rscript /tmp/packages.R
