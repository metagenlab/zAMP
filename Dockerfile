FROM continuumio/miniconda3:4.5.12

ADD envs/r_visualization.yml /tmp/r_visualization.yml
RUN conda env create -f /tmp/r_visualization.yml -n r_visualization

RUN /bin/bash -c 'source activate /opt/conda/envs/r_visualization/'

RUN apt-get update && apt-get -y install libv8-dev


RUN dpkg -L libv8-dev


RUN /bin/bash -c 'source activate /opt/conda/envs/r_visualization/ && Rscript -e "install.packages(\V8\", dependencies=TRUE, repos=\"http://cran.us.r-project.org\")"'
