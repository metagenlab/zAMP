FROM continuumio/miniconda3:4.5.12

ADD envs/r_visualization.yml /tmp/r_visualization.yml
RUN conda env create -f /tmp/r_visualization.yml -n r_visualization

RUN conda info env

# Pull the environment name out of the r_visualization.yml
#RUN echo "source activate r_visualization" > ~/.bashrc
#RUN echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc
#RUN echo "conda activate" >> ~/.bashrc
#RUN export PATH="/opt/conda/bin:$PATH"


#ENV PATH /opt/conda/envs/r_visualization/bin:$PATH

RUN /bin/bash -c 'source activate /opt/conda/envs/r_visualization/'

RUN export DEBIAN_FRONTEND=noninteractive TERM=linux && \
  apt-get update && \
  apt-get update && \
  apt-get -y --no-install-recommends install libv8-dev


RUN /bin/bash -c 'source activate /opt/conda/envs/r_visualization/ && conda info --envs'


RUN /bin/bash -c 'source activate /opt/conda/envs/r_visualization/ && Rscript -e "install.packages(randomcoloR, repos='http://cran.us.r-project.org', dependencies=TRUE)"'
