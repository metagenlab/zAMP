## Install Ubuntu
FROM ubuntu:16.04

############################## Install miniconda environement, from miniconda2 Dockerfile ##############################
#  $ docker build . -t continuumio/miniconda:latest -t continuumio/miniconda:4.5.11 -t continuumio/miniconda2:latest -t continuumio/miniconda2:4.5.11
#  $ docker run --rm -it continuumio/miniconda2:latest /bin/bash
#  $ docker push continuumio/miniconda:latest
#  $ docker push continuumio/miniconda:4.5.11
#  $ docker push continuumio/miniconda2:latest
#  $ docker push continuumio/miniconda2:4.5.11

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/conda/bin:$PATH

RUN apt-get update --fix-missing && apt-get install -y wget bzip2 ca-certificates \
    libglib2.0-0 libxext6 libsm6 libxrender1 \
    git mercurial subversion

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda2-4.5.11-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc

RUN apt-get install -y curl grep sed dpkg && \
    TINI_VERSION=`curl https://github.com/krallin/tini/releases/latest | grep -o "/v.*\"" | sed 's:^..\(.*\).$:\1:'` && \
    curl -L "https://github.com/krallin/tini/releases/download/v${TINI_VERSION}/tini_${TINI_VERSION}.deb" > tini.deb && \
    dpkg -i tini.deb && \
    rm tini.deb && \
    apt-get clean

ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD [ "/bin/bash" ]

############################## Install a default R ##############################
RUN echo "deb https://cloud.r-project.org/bin/linux/ubuntu xenial-cran35/" >> /etc/apt/sources.list && \
	apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9 && \
	apt install apt-transport-https && \
  apt update && \
	apt-get install r-base -y

## Install a dependancy for the r-V8 package, itself needed for randomcoloR
RUN apt-get update && apt-get -y install libv8-dev libcurl4-openssl-dev

############################## Import definition of the conda environment ##############################
COPY envs/r_visualization2.yml /tmp/r_visualization2.yml

############################## Activate USER ##############################
RUN useradd -r -u 1080 pipeline_user
RUN mkdir -p /home/pipeline_user
RUN chown pipeline_user -R /home/pipeline_user
USER pipeline_user

############################## Install conda env ##############################
## Create the conda environement
RUN conda env create -f /tmp/r_visualization2.yml -n r_visualization

############################## Add the needed packages ##############################
## Download the r-V8 package
RUN wget https://cran.r-project.org/src/contrib/Archive/V8/V8_1.5.tar.gz -O /tmp/rv8.tar.gz

## Install the package
RUN R CMD INSTALL /tmp/rv8.tar.gz -l /home/pipeline_user/.conda/envs/r_visualization/lib/R/library/

## Download the randomcoloR package
RUN wget https://cran.r-project.org/src/contrib/randomcoloR_1.1.0.tar.gz -O /tmp/randomcoloR.tar.gz

## Install the package
RUN R CMD INSTALL /tmp/randomcoloR.tar.gz -l /home/pipeline_user/.conda/envs/r_visualization/lib/R/library/

## Activate the r-visualiation
ENV PATH /home/pipeline_user/.conda/envs/r_visualization/bin:$PATH
WORKDIR /home/pipeline_user
RUN source activate r_visualization
