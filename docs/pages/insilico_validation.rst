
.. _insilico:

########################################################################
*In silico* validation tool
########################################################################

************************************************************************
Aim
************************************************************************

This dedicated aims at predicting *in silico* if specific taxa are:
    
    1. amplified by a set of PCR primers
    2. accurately taxonomically classified based on the generated amplicon
   

************************************************************************
Working principle
************************************************************************

Based on a user-defined list, genome assemblies are downloaded from the NCBI database with `Assembly Finder
<https://github.com/metagenlab/assembly_finder>`_. Then, PCR primers sequences provided by the user are used to run an *in silico* PCR with `Simulate_PCR <https://doi.org/10.1186/1471-2105-15-237>`. The generated in silico* amplicons are then treated by the pipeline at they would with our main pipeline (primers trimming, taxonomic classification). 

Finally, this tool provides a table with, for each of the downloaded assembly, a description of the amplicons predicted to be amplified with the tested PCR primers (number of sequence variants, number of copies, expected and obtained taxonomic classification). 


************************************************************************
Requirements
************************************************************************
This tools shares :ref:`requirements with the RSP4ABM main pipeline <_setup>` and in particular: 

- a local copy of RSP4ABM (:ref:`(cloned with git) <git>`) (with the --recursive flag to obtain Assembly Finder at the same time)
- :ref:`Snakemake <snakemake>`
- :ref:`Singularity <singularity>` (here required and not only optional)
- :ref:`a preprocessed taxonomic database <DB_preprocessing>`

In addition, the  external pipeline is required 


************************************************************************
Inputs
************************************************************************
To execute the pipeline, one needs:


1. Input table
=======================================================================

.. Note:: It is recommended to generate a new folder (outside of the pipeline itself) where these files are created and the pipeline will be executed. 

This table contains a list of taxa to be tested. First column must contain taxonomic identifiers matching the `identifiers from the NCBI taxonomy database <https://www.ncbi.nlm.nih.gov/taxonomy/>`. Alternatively, instead of their names, taxa can also be indentified by their taxID. Then, a second column must describe the number of assemblies. The two columns should be separated by a tabulation.


.. Note:: The number of assemblies from a given taxa should not exceed the number of assemblies available on the NCBI. 


**Input file example:**

.. literalinclude:: ../../ressources/template_files/16S_input_table_insilico.tsv
    :language: csv


2. Config file
=======================================================================

The config files specifies the different parameters of the pipeline as well as parameters of Assembly Finder. 

**Config file example:**

.. literalinclude:: ../../ressources/template_files/config_in_silico_validation.yaml
    :language: yaml


https://github.com/metagenlab/microbiome16S_pipeline/blob/master/ressources/template_files/16S_input_table_insilico.tsv



************************************************************************
Execution
************************************************************************

Once all the requirements installed and the input files ready, one can exectute the pipeline. In an environment where Snakemake is available, it can be run as follows: 

.. code::bash

   snakemake --snakefile {path_to_pipeline}/Insilico_taxa_assign.Snakefile  --use-singularity --singularity-prefix {path_to_singularity_images} --cores {number of cores} --configfile {path_to_config} --resources ncbi_requests={number of request to NCBI} -k

.. code::bash


