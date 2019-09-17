
### Overview
This is a comprehensive pipeline for 16S rRNA metagenomics integrating in a Snakemake worflow the best of many tools such as FASTQC, MultiQC, Qiime, rdp, EzBioCloud, Qiime2, vsearch, DADA2, phyloseq (ref to be added). It includes all steps from local .fastq or SRA depository files to generation of basic visualization plots, including quality control plots of intermediate steps. 
    
Different options are offered for use, with their pros and cons regarding reliability and flexibility. Use of a defined Snakemake pipeline and conda environments is already a good way to insure reproducibility (Method 1). However, better even reproducibility and software stability can be reached thanks to the use of containers, either with Singularity (Method 2, recommended) or Docker (Method 3).
  
### Method 1 - Snakemake with conda environnements
_Allows the more flexibility, with possibility to easily modify and personalize the pipeline. However, there are risks of errors or results inconsistency due to changes in versions. Furthermore, _simulate_PCR_ must be installed independently since it is not available through conda. 
#### Requirements:
##### Computer
A linux machine would be the best (should work as well on MacOSX, yet not tested)). At least 16Gb of RAM are needed, even more with larger datasets and depending of the used classifier. (_RDP_ requiring more RAM than _decipher_)

_Tested with Ubuntu 18.04 with 4 CPUs and 32Gb of RAM_

##### Cloned pipeline
```
git clone https://github.com/metagenlab/microbiome16S_pipeline.git
```

##### Miniconda3
Installed [following developers recommendations](https://docs.conda.io/en/latest/miniconda.html) and with relevant channels added running in a thermal the following commands :
```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set restore_free_channel true
```
_Tested with version 4.6.14_

##### Snakemake
Installed in a dedicated snakemake environments with :
```
conda create -n snakemake snakemake=5.6.0
```
_Tested with version 5.6.0_

#### Use
 ```
snakemake --snakefile ${pipeline_folder}/Snakefile --use-conda --conda-prefix ${conda_path} --cores {threads_number} --configfile {config_file_name} --resources max_copy={max_number_of_files_copied_at_the_same_time}
```
   
### Method 2 - With Singularity
_Relatively easy to set-up, and insuring software stability thanks to containerization. Here all what   dependencies are contained within the container. To be precise, we are taking advantage of the ability of Singularity to run the Docker container prepared for this pipeline. As an advantage to Docker, here no particular user settings are required_

#### Requirements:
##### Computer
As for Method 1 but here only adapted to Linux (an alpha version for Singularity exists for MacOS).

_Tested with Ubuntu 18.04 with 4 CPUs and 32Gb of RAM_

##### Singularity
Singularity is a system enabling the use of singularity or Docker containers. It should be installed [as indicated here](https://sylabs.io/guides/3.1/user-guide/quick_start.html#quick-installation-steps). Alternatively, it can be [installed through Conda](https://anaconda.org/conda-forge/singularity). However, in this last case, one need to adapt the file ownership and access rights as indicated [here](https://github.com/conda-forge/singularity-feedstock/blob/master/recipe/post-link.sh).  
  
  _Tested with version 3.0.1_

#### Commands:

```
## Run the contained interactively
singularity run -e -H /home/pipeline_user/ docker://metagenlab/amplicon_pipeline:v.0.9.9

## Run the pipeline from into the cointainer.
snakemake --snakefile /home/pipeline_user/microbiome16S_pipeline/Snakefile --use-conda --conda-prefix /opt/conda/ --cores {threads_number} --configfile {config_file_path} --resources max_copy={max_number_of_files_copied_at_the_same_time} mem_mb = {available_memory}
```


### Method 3 - With Docker
##### Computer
Works on Windows, MacOS and Linux.
_Tested on Linux, Windows 10 and MacOSX_

##### User settings: 
Our Docker image is fitted for a user called "pipeline_user" whose UID is 1080. It is advised to create this user on your computer before using the Docker image to run your analysis:

```
sudo useradd -G docker,sudo -u 1080 pipeline_user
sudo mkdir /home/pipeline_user/
sudo chown pipeline_user -R /home/pipeline_user/
sudo passwd pipeline_user
```

Alternatively, you can run the Docker as root (--user root) but the created folders will belong to the root user of your computer.

##### Docker
Install the CE version following these [instructions](https://docs.docker.com/install/linux/docker-ce/ubuntu/) for ubuntu. Also make sure you have created the docker group and that you can run Docker without sudo following these [instruction](https://docs.docker.com/install/linux/linux-postinstall/). If you can't have access to the internet when inside a Docker container, apply those [changes](https://docs.docker.com/install/linux/linux-postinstall/#disable-dnsmasq).
  
  
#### Use
_Connected as pipeline_user :_ 
```
docker run -it --rm --mount source="$(pwd)",target=/home/pipeline_user/data/analysis/,type=bind metagenlab/amplicon_pipeline:v.0.9.9
```
and then
```
snakemake --snakefile /home/pipeline_user/microbiome16S_pipeline/Snakefile --use-conda --conda-prefix /opt/conda/ --cores {threads_number} --configfile {config_file_path} --resources max_copy={max_number_of_files_copied_at_the_same_time}
```
or directly

```
docker run -it --rm --mount source="$(pwd)",target=/home/pipeline_user/data/analysis/,type=bind metagenlab/amplicon_pipeline:v.0.9.9 \ sh -c 'snakemake --snakefile /home/pipeline_user/microbiome16S_pipeline/Snakefile --use-conda --conda-prefix /opt/conda/ --cores {threads_number} --configfile {config_file_path} --resources max_copy={max_number_of_files_copied_at_the_same_time}
```



