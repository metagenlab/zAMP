
### Overview
This is a comprehensive pipeline for amplicon-based metagenomics integrating in a Snakemake workflow the [best functions of many tools](#References). It enables performant and reproducibile processing of 16S rRNA or ITS Illumina paired-end reads. The whole process from local .fastq or SRA depository files to generation of basic visualization plots, including quality control plots of intermediate steps, is covered. 
    
The Snakemake pipeline can be exectued by cloning this repository and relying on conda environments (Method 1). For better reproducibility, software stability and convenient installation, our stand-alone container can be used, either with Singularity (Method 2, recommended) or Docker (Method 3).
  
### Method 1 - Snakemake with conda environnements
_Allows flexibility, with possibility to easily modify and personalize the pipeline. However, there are risks of errors or result inconsistencies due to changes in versions. Furthermore, _simulate_PCR_ must be installed independently for the in silico validation, since it is not available through conda._ 
#### Requirements:
##### Computer
A linux machine would be the best (should work as well on MacOSX, yet not tested). At least 16Gb of RAM are needed, even more with larger datasets and depending of the used classifier. (_RDP_ requiring more RAM than _decipher_)

_Tested with Ubuntu 18.04 with 4 CPUs and 32Gb of RAM_

##### Cloned pipeline
```
git clone https://github.com/metagenlab/microbiome16S_pipeline.git
```

##### Miniconda3
Installed [following developers' recommendations](https://docs.conda.io/en/latest/miniconda.html) and with relevant channels added running in a thermal the following commands :
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
_Relatively easy to set-up and easier to work with than Docker, due to simpler user-rights management. We take advantage of the ability of Singularity to run the Docker container prepared for this pipeline. Insures software stability thanks to containerization. Here, all dependencies are contained within the container._

#### Requirements:
##### Computer
As for Method 1 but here only adapted to Linux (an alpha version for Singularity exists for MacOS).

_Tested with Ubuntu 18.04 with 4 CPUs and 32Gb of RAM_

##### Singularity
Singularity is a system enabling the use of singularity or Docker containers. It should be installed [as indicated here](https://sylabs.io/guides/3.1/user-guide/quick_start.html#quick-installation-steps). 
  
  _Tested with version 3.0.1_

#### Commands:

```

## Run the container interactively
singularity shell docker://metagenlab/amplicon_pipeline:v.0.9.10

## Run the pipeline from within the container.
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
docker run -it --rm --mount source="$(pwd)",target=/home/pipeline_user/data/analysis/,type=bind metagenlab/amplicon_pipeline:v.0.9.10
```
and then
```
snakemake --snakefile /home/pipeline_user/microbiome16S_pipeline/Snakefile --use-conda --conda-prefix /opt/conda/ --cores {threads_number} --configfile {config_file_path} --resources max_copy={max_number_of_files_copied_at_the_same_time}
```
or directly

```
docker run -it --rm --mount source="$(pwd)",target=/home/pipeline_user/data/analysis/,type=bind metagenlab/amplicon_pipeline:v.0.9.10 \ sh -c 'snakemake --snakefile /home/pipeline_user/microbiome16S_pipeline/Snakefile --use-conda --conda-prefix /opt/conda/ --cores {threads_number} --configfile {config_file_path} --resources max_copy={max_number_of_files_copied_at_the_same_time}
```



### References 
- Snakemake

    *Köster, J., & Rahmann, S. (2012). Snakemake-a scalable bioinformatics workflow engine. Bioinformatics, 28(19), 2520–2522. https://doi.org/10.1093/bioinformatics/bts480*

- FASTQC

    *Andrews, S. (2010). FASTQC. A quality control tool for high throughput sequence data. 2010. Http://Www.Bioinformatics.Babraham.Ac.Uk/Projects/Fastqc/*

- MultiQC

    *Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics, 32(19), 3047–3048. https://doi.org/10.1093/bioinformatics/btw354*

- DADA2

    *Callahan, B. J., McMurdie, P. J., Rosen, M. J., Han, A. W., Johnson, A. J. A., & Holmes, S. P. (2016). DADA2: High-resolution sample inference from Illumina amplicon data. Nature Methods, 13(7), 581–583. https://doi.org/10.1038/nmeth.3869*

- VSEARCH

    *Rognes, T., Flouri, T., Nichols, B., Quince, C., & Mahé, F. (2016). VSEARCH: a versatile open source tool for metagenomics. PeerJ, 4, e2584. https://doi.org/10.7717/peerj.2584*


- Qiime

    *Caporaso, J. G., Kuczynski, J., Stombaugh, J., Bittinger, K., Bushman, F. D., Costello, E. K., … Knight, R. (2010). QIIME allows analysis of high-throughput community sequencing data. Nature Methods, 7(5), 335–336. https://doi.org/10.1038/nmeth.f.303*

- Qiime2 

    *Bolyen, E., Dillon, M., Bokulich, N., Abnet, C., Al-Ghalith, G., Alexander, H., … Caporaso, G. (2018). QIIME 2: Reproducible, interactive, scalable, and extensible microbiome data science. PeerJ Preprints. https://doi.org/10.7287/peerj.preprints.27295*

- RDP

    *Wang, Q., Garrity, G. M., Tiedje, J. M., & Cole, J. R. (2007). Naive Bayesian classifier for rapid assignment of rRNA sequences into the new bacterial taxonomy. Applied and Environmental Microbiology, 73(16), 5261–5267. https://doi.org/10.1128/AEM.00062-07*
    
- IDTAXA in Decipher

    *Murali, A., Bhargava, A., & Wright, E. S. (2018). IDTAXA: A novel approach for accurate taxonomic classification of microbiome sequences. Microbiome. https://doi.org/10.1186/s40168-018-0521-5*

- EzBioCloud

    *Yoon, S.-H., Ha, S.-M., Kwon, S., Lim, J., Kim, Y., Seo, H., & Chun, J. (2017). Introducing EzBioCloud: a taxonomically united database of 16S rRNA gene sequences and whole-genome assemblies. International Journal of Systematic and Evolutionary Microbiology, 67(5), 1613–1617. https://doi.org/10.1099/ijsem.0.001755*

- Silva 

    *Quast, C., Pruesse, E., Yilmaz, P., Gerken, J., Schweer, T., Yarza, P., … Glöckner, F. O. (2013). The SILVA ribosomal RNA gene database project: Improved data processing and web-based tools. Nucleic Acids Research. https://doi.org/10.1093/nar/gks1219*

- phyloseq

    *McMurdie, P. J., & Holmes, S. (2013). phyloseq: An R Package for Reproducible Interactive Analysis and Graphics of Microbiome Census Data. PLoS ONE, 8(4), e61217. https://doi.org/10.1371/journal.pone.0061217*

- Krona 

    *Ondov, B. D., Bergman, N. H., & Phillippy, A. M. (2011). Interactive metagenomic visualization in a Web browser. BMC Bioinformatics, 12(1), 385. https://doi.org/10.1186/1471-2105-12-385*

- ALDex2

    *Fernandes, A. D., Reid, J. N. S., Macklaim, J. M., McMurrough, T. A., Edgell, D. R., & Gloor, G. B. (2014). Unifying the analysis of high-throughput sequencing datasets: Characterizing RNA-seq, 16S rRNA gene sequencing and selective growth experiments by compositional data analysis. Microbiome, 2(1), 1–13. https://doi.org/10.1186/2049-2618-2-15*

- Vegan
    
    *Oksanen, J., Kindt, R., Legendre, P., O’Hara, B., Simpson, G. L., Solymos, P. M., … & Wagner, H. (2008). The vegan package. Community Ecology Package, (May 2014), 190. Retrieved from https://bcrc.bio.umass.edu/biometry/images/8/85/Vegan.pdf*

- metagenomeSeq

    *Joseph, A., Paulson, N., Olson, N. D., Wagner, J., Talukder, H., & Corrada, H. (2019). Package ‘ metagenomeSeq .’*
    
- edgeR

    *Robinson, M. D., McCarthy, D. J., & Smyth, G. K. (2009). edgeR: A Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics, 26(1), 139–140. https://doi.org/10.1093/bioinformatics/btp616*
    
- Simulate_PCR

    *Gardner, S. N., & Slezak, T. (2014). Simulate_PCR for amplicon prediction and annotation from multiplex, degenerate primers and probes. BMC Bioinformatics, 15(1), 2–7. https://doi.org/10.1186/1471-2105-15-237*

