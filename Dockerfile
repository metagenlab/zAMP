## Install Ubuntu
FROM ubuntu:18.04

############################## Install miniconda environement, from miniconda3 Dockerfile ##############################
#  $ docker build . -t continuumio/miniconda3:latest -t continuumio/miniconda3:4.5.11
#  $ docker run --rm -it continuumio/miniconda3:latest /bin/bash
#  $ docker push continuumio/miniconda3:latest
#  $ docker push continuumio/miniconda3:4.5.11

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/conda/bin:$PATH
ENV TZ Europe/Zurich

RUN echo $TZ > /etc/timezone && \
    apt-get update && apt-get install -y tzdata && \
    rm /etc/localtime && \
    ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && \
    dpkg-reconfigure -f noninteractive tzdata && \
    apt-get clean

RUN apt-get update --fix-missing && \
    apt-get install -y wget nano bzip2 ca-certificates curl git && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-4.6.14-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    /opt/conda/bin/conda clean -tipsy && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc

ENV TINI_VERSION v0.16.1
ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /usr/bin/tini
RUN chmod +x /usr/bin/tini

############################## Create pipeline_user, set useful variables ##############################
RUN useradd -r -u 1080 pipeline_user
ENV main=/home/pipeline_user
WORKDIR $main
ENV pipeline_folder=${main}/microbiome16S_pipeline

########################### Install java (needed for Qiime tax assignemnt) and Snakemake ##############################
RUN conda config --add channels defaults && conda config --add channels bioconda && conda config --add channels conda-forge
RUN conda install snakemake=5.5.4 java-jdk perl-bioperl conda=4.6.14

######### Install PANDAseq (libltdl7) and r-v8 (libv8-dev) dependances, and a package required for png plotting  (libcairo2) ##########
RUN apt-get update && apt-get install libltdl7 libv8-dev libcairo2-dev -y

############################## Get the pipeline through github #######################
## Call the access token to reach the private github repo
ARG GITHUB_AT

## Clone the github
RUN git clone --single-branch --branch master https://$GITHUB_AT@github.com/metagenlab/microbiome16S_pipeline.git $pipeline_folder

## cd in the validation directory
WORKDIR ${pipeline_folder}/data/validation_datasets

#################### Build environements of the pipeline #####################
## Here, with "--create-envs-only", we only build the environements
RUN snakemake --snakefile ${pipeline_folder}/Snakefile --use-conda --conda-prefix /opt/conda/ --create-envs-only --configfile config.yml all PICRUSt2_output

##################### Install r-v8 dependancy, r-v8 and randomcoloR R package, used for plotting ######################
### Warning:
#### - we need "libgcc-ng=7.2.0" installed in the environement for the R-V8 installation to work
#### - libv8-dev must be installed and the containing path (--configure-vars=\"INCLUDE_DIR=/usr/include/ LIB_DIR=/usr/lib/ \") correctly pointed at
#### - Due to conda dependacies definition and incombatility between "libgcc-ng = 7.2.0" needed by R-V8 to build and this build, we cannot the "build=*_1007" in the env definition. However, this build prevent "cannot access ldpath error". Hences, we install it after the facts.
## Download packages
#RUN wget https://cran.r-project.org/src/contrib/V8_2.3.tar.gz -O /tmp/rv8.tar.gz
RUN wget https://cran.r-project.org/src/contrib/randomcoloR_1.1.0.tar.gz -O /tmp/randomcoloR.tar.gz

## Recover the specific env by its name and update the env with randomcoloR
#&& R CMD INSTALL --configure-vars=\"INCLUDE_DIR=/usr/include/ LIB_DIR=/usr/lib/ \" /tmp/rv8.tar.gz &&
RUN export barplots_env=$(basename $(grep "name: barplots" /opt/conda/*.yaml | cut -d: -f1) .yaml) && /bin/bash -ce """source activate /opt/conda/${barplots_env} && Rscript -e \"install.packages('randomcoloR', repos = 'https://stat.ethz.ch/CRAN/')\" && conda install r-base==3.5.1[build=*_1007] """
RUN rm /tmp/randomcoloR.tar.gz

## Install simulate PCR, DOI: 10.1186/1471-2105-15-237 for amplicons validation
RUN wget --quiet https://sourceforge.net/projects/simulatepcr/files/simulate_PCR-v1.2.tar.gz/download -O simulate_PCR.tar.gz && mkdir /opt/simulate_PCR && tar xzf simulate_PCR.tar.gz -C /opt/simulate_PCR && rm simulate_PCR.tar.gz
ENV PATH="/opt/simulate_PCR:${PATH}"
RUN apt-get install libgd-dev &&  cpanm GD && cpanm 'XML::DOM::XPath' --force && cpanm Bio::SeqIO

################# Clean unnecessary packages ###################
RUN conda clean -a
RUN apt-get autoremove -y

#################### Run the pipeline to test it #####################
## Here, we run the pipeline to test it, without PICRUST as output since it is computationally very demanding
RUN snakemake --snakefile ${pipeline_folder}/Snakefile --cores 50 --resources max_copy=1 --use-conda --conda-prefix /opt/conda/ --configfile ${pipeline_folder}/data/validation_datasets/config.yml all

#################### Set final access rights and working dir #####################
RUN chown -R pipeline_user ${main}/
USER pipeline_user
RUN mkdir -p ${main}/data/analysis/
RUN conda init bash
WORKDIR ${main}/data/analysis/
ENTRYPOINT ["/bin/bash"]
