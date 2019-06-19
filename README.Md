
### Overview
This is a comprehensive pipeline for 16S rRNA metagenomics integrating the best of many tools such as FASTQC, MultiQC, Qiime, rdp, EzBioCloud, Qiime2, vsearch, DADA2, phyloseq (ref to be added) in a Snakemake worflow. It includes all steps from .fastq recovery through symbolic links or SRA depository to generation of basic visualization plots, including quality control plots at different stages. 
    
  There are various ways to use the pipeline. The use of a defined Snakemake pipeline already enhance reproducibility (method 1). However, best reproducibility is reached when using Singularity (method 2) or Docker (method 3) container-based approach.
  
### Method 1 - Snakemake in conda environnements
_Allows the best flexibility with modification of the pipeline_

#### Requirements:
##### Computer
A linux machine would be the best,  however, should work as well on MacOSX (to be tested). At least 16Gb of RAM are needed, even more with larger datasets.

_Tested with Ubuntu 18.04 with 4 CPUs and 32Gb of RAM_

#### Cloned pipeline
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
Installed by :
```
conda create -n snakemake snakemake=5.5.1
```


_Tested with version 5.5.1_

##### Singularity (optional)
_Singularity is a system enabling the use of docker containers for per-rule environment [see Snakemake doc](https://snakemake.readthedocs.io/en/v5.5.1/snakefiles/deployment.html). Library needed by the "barplots" rule cannot be installed through conda. As a workaround, we have built a Docker image which is runned by singularity when the "--use-singularity" flag is present in the snakemake command._
   
  
  It must be installed [as indicated here](https://www.sylabs.io/guides/3.0/user-guide.pdf). 
  
  _Tested with version 3.0.1_
  
#### Use:
##### With "Barplots: False" in config (no singularity):
```
snakemake --snakefile ${pipeline_folder}/Snakefile --use-conda --conda-prefix ${conda_path} --cores {threads_number} --configfile {config_file_path} --resources max_copy={max_number_of_files_copied_at_the_same_time}
```

##### With "Barplots: True" in config (incl. singularity):
 ```
snakemake --snakefile ${pipeline_folder}/Snakefile --use-conda --conda-prefix ${conda_path} --cores {threads_number} --configfile {config_file_name} --resources max_copy={max_number_of_files_copied_at_the_same_time} --use-singularity
```
 
  
### Method 2 - With Singularity
_Relatively easy to set-up and insuring reproducibility thanks to containerization. Here all what is needed to run the pipeline is a container. In fact, we are taking advantage of the abilify of Singularity to run the Docker container prepared for this pipeline, with the added benefits of solving the user access problems (see Method 3) and allowing the access to all local files._

#### Requirements:
##### Computer
As for Method 1 but here normally only for linux. However, an alpha version for Singularity exists for MacOS.

##### Singularity
As for Method 1

#### Use:

```
## Run the contained interactively
singularity run docker://valscherz/amplicon_snakemake_pipeline:0.0.9-beta

## Run the pipeline
snakemake --snakefile /home/pipeline_user/microbiome16S_pipeline/Snakefile --use-conda --conda-prefix /opt/conda/ --cores {threads_number} --configfile {config_file_path} --resources max_copy={max_number_of_files_copied_at_the_same_time}
```
  

### Method 3 - With Docker
##### Computer
Should work on Windows, MacOS and Linux.
_Tested on Linux_

##### User settings: 
Our Docker image is fit for a user called pipeline_user whose UID is 1080. It is advised to create this user on your computer before using the Docker image to run your analysis:

```
sudo useradd -G docker,sudo -u 1080 pipeline_user
sudo mkdir /home/pipeline_user/
sudo chown pipeline_user -R /home/pipeline_user/
sudo passwd pipeline_user
```

Alternatively, you can run the Docker as root (--user root) but the created folders will belong to the root user of your computer.

##### Docker
Install the CE version following these [instructions](https://docs.docker.com/install/linux/docker-ce/ubuntu/) for ubuntu. Also make sure you have created the docker group and that you can run docker without sudo following these [instruction](https://docs.docker.com/install/linux/linux-postinstall/). If you can't have access to the internet when inside a Docker container, apply those [changes](https://docs.docker.com/install/linux/linux-postinstall/#disable-dnsmasq).
  
  
## General use
```
snakemake --snakefile $pipeline_folder/Snakefile --use-conda --conda-prefix $pipeline_folder/miniconda3/ --cores 8 --configfile config.yaml --use-singularity -k <disered_output_rule>
```


_where:_

_--use-singularity_ to allow the use of singularity, a system enabling the use of docker environments for certain rules

_--k_ to go on with the pipeline despite errors

