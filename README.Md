
### Overview
This is a comprehensive pipeline for 16S rRNA metagenomics integrating in a Snakemake worflow the best of many tools such as FASTQC, MultiQC, Qiime, rdp, EzBioCloud, Qiime2, vsearch, DADA2, phyloseq (ref to be added). It includes all steps from local .fastq or SRA depository to generation of basic visualization plots, including quality control plots of intermediate steps. 
    
There are various ways to use the pipeline, with their pros and cons regarding fiability and flexbility. Use of a defined Snakemake pipeline and conda environments is already a good way to insure reproducibility (Method 1). However, better even reproducibility and avoidance of softwares compatibility issues can be reached thanks to the use of containers, either with Singularity (Method 2, recommended) or Docker (Method 3).
  
### Method 1 - Snakemake with conda environnements
_Allows the best flexibility, with possibility to easily modify and personalize the pipeline. However, there are risks of errors or results inconsistency due to changes in versions. Furthermore, the "barplots" option in config must be set on "False", or the randomcolorR package muste installed manually and made avaible in the "barplots" environements, since this package in required by the Barplots rule but is not available in conda_

#### Requirements:
##### Computer
A linux machine would be the best, but should work as well on MacOSX (to be tested, maybe conda environement definitions to be adated). At least 16Gb of RAM are needed, even more with larger datasets.

_Tested with Ubuntu 18.04 with 4 CPUs and 32Gb of RAM_

##### Cloned pipeline
```
git clone https://github.com/metagenlab/microbiome16S_pipeline.git
```

##### Miniconda3
Installed [following developpers recommendations](https://docs.conda.io/en/latest/miniconda.html) and with relevent channels added runnning in a termal the following commands :
```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```
_Tested with version 4.6.14_

##### Snakemake
Installed in a dedicated snakemake environements with :
```
conda create -n snakemake snakemake=5.5.1
```
_Tested with version 5.5.1_

#### Use
 ```
snakemake --snakefile ${pipeline_folder}/Snakefile --use-conda --conda-prefix ${conda_path} --cores {threads_number} --configfile {config_file_name} --resources max_copy={max_number_of_files_copied_at_the_same_time}
```
   
### Method 2 - With Singularity
_Relatively easy to set-up and insuring reproducibility thanks to containerization. Here all what is needed to run the pipeline is Singularity installed. To be precise, we are taking advantage of the abilify of Singularity to run the Docker container prepared for this pipeline. The advantage compared with Docker is that with Singularity not specific user is required._

#### Requirements:
##### Computer
As for Method 1 but here normally only for linux. However, an alpha version for Singularity exists for MacOS.

_Tested with Ubuntu 18.04 with 4 CPUs and 32Gb of RAM_

##### Singularity
Singularity is a system enabling the use of singularity or docker containers. It should be installed [as indicated here](https://www.sylabs.io/guides/3.0/user-guide.pdf). (Installation through apt or conda could also work be should be tested.)
  
  _Tested with version 3.0.1_

#### Commands:

```
## Run the contained interactively
singularity run docker://metagenlab/amplicon_pipeline:v.0.9.6-beta

## Run the pipeline from into the cointainer.
snakemake --snakefile /home/pipeline_user/microbiome16S_pipeline/Snakefile --use-conda --conda-prefix /opt/conda/ --cores {threads_number} --configfile {config_file_path} --resources max_copy={max_number_of_files_copied_at_the_same_time} mem_mb = {available_memory}
```


### Method 3 - With Docker
##### Computer
Should work on Windows, MacOS and Linux.
_Tested on Linux_

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
Install the CE version following these [instructions](https://docs.docker.com/install/linux/docker-ce/ubuntu/) for ubuntu. Also make sure you have created the docker group and that you can run docker without sudo following these [instruction](https://docs.docker.com/install/linux/linux-postinstall/). If you can't have access to the internet when inside a Docker container, apply those [changes](https://docs.docker.com/install/linux/linux-postinstall/#disable-dnsmasq).
  
  
#### Use
```
docker run -it --rm --mount source="$(pwd)",target=/home/pipeline_user/data/analysis/,type=bind metagenlab/amplicon_pipeline:v.0.9.6-beta
```
and then
```
snakemake --snakefile /home/pipeline_user/microbiome16S_pipeline/Snakefile --use-conda --conda-prefix /opt/conda/ --cores {threads_number} --configfile {config_file_path} --resources max_copy={max_number_of_files_copied_at_the_same_time}
```
or directly

```
docker run -it --rm --mount source="$(pwd)",target=/home/pipeline_user/data/analysis/,type=bind metagenlab/amplicon_pipeline:v.0.9.6-beta \ sh -c 'snakemake --snakefile /home/pipeline_user/microbiome16S_pipeline/Snakefile --use-conda --conda-prefix /opt/conda/ --cores {threads_number} --configfile {config_file_path} --resources max_copy={max_number_of_files_copied_at_the_same_time}
```



