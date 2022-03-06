
.. _insilico:

########################################################################
*In silico* validation tool
########################################################################

************************************************************************
Aim
************************************************************************

This dedicated aims at predicting *in silico* if specific taxa are:
    
    1. amplified by a set of PCR primers used for amplicon-based metagenomics
    2. accurately taxonomically classified based on the generated amplicon
   

************************************************************************
Working principle
************************************************************************

Based on a user-defined list, genome assemblies are downloaded from the NCBI database with `Assembly Finder <https://github.com/metagenlab/assembly_finder>`_. Then, PCR primer sequences provided by the user are used to run an *in silico* PCR with `Simulate_PCR <https://doi.org/10.1186/1471-2105-15-237>`_. To the their discriminative power, the generated *in silico* amplicons are  treated by the pipeline as they would if they were the results of sequencing reads (primers trimming, taxonomic classification). 

Finally, this tool provides a table with, for each of the downloaded assembly, a description of the amplicons predicted to be amplified with the PCR primers (number of sequence variants, number of copies, expected and obtained taxonomic classification). 


************************************************************************
Requirements
************************************************************************
This tools shares :ref:`requirements with the RSP4ABM main pipeline <_setup>`. Thus, it requires: 

- a local copy of RSP4ABM (:ref:`(cloned with git) <git>`) (with the --recursive flag to obtain Assembly Finder at the same time)
- :ref:`Snakemake <snakemake>`
- :ref:`Singularity <singularity>` (here required and not only optional)
- :ref:`A taxonomic database preprocessed with our dedicated pipeline <DB_preprocessing>`


************************************************************************
Inputs
************************************************************************
To execute the pipeline, one needs:


1. Input table
=======================================================================

.. Note:: It is recommended to generate a new folder (outside of the pipeline itself) where these files are created and the pipeline executed. 

This table contains a list of taxa to be tested. First column must contain taxonomic identifiers matching the `identifiers from the NCBI taxonomy database <https://www.ncbi.nlm.nih.gov/taxonomy/>`_. Alternatively, instead of their names, taxa can also be indentified by their taxID. Then, a second column must describe the number of assemblies. The two columns should be separated by a tabulation.


.. Hint:: The number of assemblies from a given taxa should not exceed the number of assemblies available on the NCBI. 


**Input file example:**

.. literalinclude:: ../../ressources/template_files/16S_input_table_insilico.tsv
    :language: csv


2. Config file
=======================================================================

The config files specifies the different parameters of the pipeline as well as parameters for `Assembly Finder <https://github.com/metagenlab/assembly_finder>`_. 

**Config file example:**

.. literalinclude:: ../../ressources/template_files/config_in_silico_validation.yaml
    :language: yaml


https://github.com/metagenlab/microbiome16S_pipeline/blob/master/ressources/template_files/16S_input_table_insilico.tsv



************************************************************************
Execution
************************************************************************

Once all the requirements installed and the input files ready, one can exectute the pipeline. In an environment where :ref:`Snakemake <snakemake>` is available, it can be run as follows: 

.. code-block:: console

    snakemake --snakefile {path_to_pipeline}/Insilico_taxa_assign.Snakefile  --use-singularity --singularity-prefix {path_to_singularity_images} --cores {number of cores} --configfile {path_to_config} --resources ncbi_requests={number of request to NCBI} -k


 
