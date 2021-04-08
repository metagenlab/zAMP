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

- a local copy of RSP4ABM (:ref:`(cloned with git) <git>`)
- :ref:`Snakemake <snakemake>`
- :ref:`Singularity <singularity>` (here required and not only optional)
- :ref:`a preprocessed taxonomic database <DB_preprocessing>`

In addition, the  external pipeline is required 





Title3
=======================================================================

Title4
-----------------------------------------------------------------------

