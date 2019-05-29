## Requirements:
#### miniconda3
_With export PATH="<path/to>/miniconda3/bin:$PATH" in .bashrc. *Careful: for some reason, taxonomic assignment with Qiime can fail is this pattern is not respected.*_

#### snakemake
_version >5.3.0 installed and containing environment activated_

#### Singularity
_Singularity is a system enabling the use of docker containers for per-rule environment [see Snakemake doc](https://snakemake.readthedocs.io/en/v5.4.0/snakefiles/deployment.html) when plotting is desired. It must be installed [as indicated here](https://www.sylabs.io/guides/3.0/user-guide.pdf). It is used here since R randomcoloR package, needed for plotting, is not available for Conda due to a missing dependancy, r-V8. ([See here for more about failure in building the package](https://github.com/conda-forge/staged-recipes/pull/6982]))._

## General use
```
snakemake --snakefile $pipeline_folder/Snakefile --use-conda --conda-prefix $pipeline_folder/miniconda3/ --cores 8 --configfile config.yaml --use-singularity -k <disered_output_rule>
```


_where:_

_--use-singularity_ to allow the use of singularity, a system enabling the use of docker environments for certain rules

_--k_ to go on with the pipeline despite errors

