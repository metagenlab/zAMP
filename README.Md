## Requirements:
- miniconda3**
** With export PATH="<path/to>/miniconda3/bin:$PATH" in .bashrc. * Careful: for some reason, it taxonomic assignment with Qiime fails if "" are missing * **
- snakemake
** version >5.3.0 installed and containing environement activaited **

## General use
```
snakemake --snakefile $pipeline_folder/Snakefile --use-conda --conda-prefix $pipeline_folder/miniconda3/ --cores 8 --configfile config.yaml --use-singularity -k <disered_output>
```


**where:

/--use-singularity/ to allow the use of singularity, a system enabling the use of docker environements for certain rules

/--k/ to go on with the pipeline despite errors

*Singularity:*

Singularity is a system enabling the use of docker containers for per-rule environements [[https://snakemake.readthedocs.io/en/v5.4.0/snakefiles/deployment.html][(see Snakemake doc)]]. It must be installed [[https://www.sylabs.io/guides/3.0/user-guide.pdf][as indicated here]]. It is used here since R randomcoloR package, needed for plotting, is not available for Conda due to a missing dependancy, r-V8. ([[https://github.com/conda-forge/staged-recipes/pull/6982][See here for more about failure in building the package]]).
