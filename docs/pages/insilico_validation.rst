
.. _insilico:

########################################################################
*In silico* validation tool
########################################################################

************************************************************************
Aim
************************************************************************

This module aims at predicting *in silico* if specific taxa are:
    
    1. amplified by a set of PCR primers used for amplicon-based metagenomics
    2. accurately classified taxonomically based on the generated amplicon
   

************************************************************************
Working principle
************************************************************************

Based on a user-defined list of tax IDs or species names, genome assemblies are downloaded from the NCBI database with `Assembly Finder <https://github.com/metagenlab/assembly_finder>`_. Then, PCR primer sequences provided by the user are used to run an *in silico* PCR with `dicey <https://github.com/gear-genomics/dicey/>`. The generated *in silico* amplicons are  treated by the pipeline as they would if they were the results of sequencing reads (primer trimming, amplicon clustering into representative sequences, taxonomic classification). 

Finally, for each of the downloaded assembly, this module provides a table with a description of the amplicons predicted to be amplified with the PCR primers (number of sequence variants, number of copies) as well as the expected and obtained taxonomic assignment. 


************************************************************************
Requirements
************************************************************************
This tools shares :ref:`requirements with the zAMP main pipeline <_setup>`. Thus, it requires: 

- a local copy of zAMP (:ref:`(cloned with git) <git>`) (with the --recursive flag to obtain Assembly Finder at the same time)
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

This table contains a list of taxa to be tested. This is a one-column table named "UserInputNames" and contains taxonomic identifiers (taxID) matching the `identifiers from the NCBI taxonomy database <https://www.ncbi.nlm.nih.gov/taxonomy/>`_. Alternatively, taxa can also be indentified by their names.  

**Input file example:**

.. literalinclude:: ../../ressources/template_files/insilico_ITS_input_taxID.tsv
    :language: csv


2. Config file
=======================================================================

The config files specifies the different parameters of the pipeline as well as parameters for `Assembly Finder <https://github.com/metagenlab/assembly_finder>`_. 

**Config file example:**

.. literalinclude:: ../../ressources/template_files/config_in_silico_validation.yaml
    :language: yaml



************************************************************************
Execution
************************************************************************

Once all the requirements installed and the input files ready, one can exectute the pipeline. In an environment where :ref:`Snakemake <snakemake>` is available, it can be run as follows: 

.. code-block:: console

    snakemake --snakefile {path_to_pipeline}/Insilico_taxa_assign.Snakefile  --use-singularity --singularity-prefix {path_to_singularity_images} --cores {number of cores} --configfile {path_to_config} --resources ncbi_requests={number of request to NCBI} -k


Alternatively, using conda:

.. code-block:: console

    snakemake --snakefile {path_to_pipeline}/Insilico_taxa_assign.Snakefile  --use-conda --cores {number of cores} --configfile {path_to_config} --resources ncbi_requests={number of request to NCBI} -k


************************************************************************
Output
************************************************************************
The pipeline gathers information on available assemblies for the requested taxIDs in "tables" folder, and the downloaded assemblies in "assemblies_gz".

The output of the in-silico amplification is in Insilico folder, and contains the following subfolders:

- PCR: contains the output of in-silico PCR amplification
- 1a_trimmed_primers: contains output of primer trimming with cutadapt
- 1c_derep: all amplified products combined into one fasta file (similar to the input that would be used in the main pipeline)
- 2_denoised: output of clustering and denoising into representative sequences, and count tables
- 3_classified: output of taxonomic classification and tables comparing expected and obtained taxonomic assignations (InSilico_compare_tax.tsv and InSilico_compare_tax_long.tsv)
