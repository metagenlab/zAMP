## DB processing :

Preparation of reference databases for our amplicon-based metagenomic snakemake pipeline is here embedded as an ad-hoc snakemake pipeline. It takes as input a reference database containing a fasta file (sequences) and a taxonomy file (corresponding taxonomy) in the Qiime format. 

**This ad-hoc reference database preparation pipeline relies on  :**
 - A dedicated snakefile (Snakefile_DB) to run this ad-hoc pipeline.
 - Configuration file (config_DB.yaml) to point at the input and output paths and defines parameter for amplicons extraction and taxonomy formatting. 

In brief, through the few rules in "EzBioCloud.rules", we will import the taxonomy and sequences in Qiime2 ".qza" format, and extract the region targeted by the PCR primers defined in the config file using Qiime2 "qiime feature-classifier extract-reads" function. Then, the sequences will be dereplicated (to keep only unique sequences) using vsearch. Following dereplication, taxonomy description will be merged to describe the different taxa that could correspond to a given sequences. A parameter in config allows to define up to how many Genus/Species names will merged. If there are more taxa matching a sequence, then a spaceholder (".sp") will be used. The number of taxa matching a specific taxonomic ID will also be written. 


**To compile a new database :**
- Create a new directory in "/data" using "{dbname}{version}.{preparationyearmonth}
- In this directory, include a Readme.Md file describing the origin of the reference database.
- Locate the fasta, taxonomy files and the desired output path (just created directory) in configfile. 
- Run the following command from the just created directory.
``` 
snakemake --snakefile $pipeline_folder/Snakefile_DB --use-conda --conda-prefix $condaprefix --cores 1  --configfile config_DB.yaml     
```
- Add to git the output fasta and taxonomy as well as the configfile, the generated logs and the Readme describing the origin of the database to enable traceability. 
