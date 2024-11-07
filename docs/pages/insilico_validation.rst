
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

Based on a user-defined list of NCBI tax IDs, assemblies or taxon queries, genome assemblies are downloaded from the NCBI database with `Assembly Finder <https://github.com/metagenlab/assembly_finder>`_. Then, PCR primer sequences provided by the user are used to run an *in silico* PCR with `in_silico_pcr <https://github.com/egonozer/in_silico_pcr>`_ (or alternatively, with `simulate_PCR <https://github.com/metagenlab/updated_simulate_PCR>`_. The generated *in silico* amplicons are  treated by the main pipeline as they would if they were the results of sequencing reads (primer trimming, amplicon clustering into representative sequences, taxonomic classification). 

Finally, for each of the downloaded assembly, this module provides a table with a description of the amplicons predicted to be amplified with the PCR primers (number of sequence variants, number of copies) as well as the expected and obtained taxonomic assignment. 


************************************************************************
Inputs
************************************************************************
To execute the pipeline, one needs:

* An input file containing the accession names or the Tax IDs of interest. This is a one-column text file without headers. The identifiers should match NCBI taxonomy. One can skip this text file and use a query term instead, see usage cases below.

* :ref:`A taxonomic database preprocessed with our dedicated pipeline <DB_preprocessing>`


**Input file example:**

.. literalinclude:: ../zamp/data/bacteria-accs.txt
    :language: csv

.. literalinclude:: ../zamp/data/fungi-taxa.txt
    :language: csv

************************************************************************
Execution
************************************************************************
Example usage cases:

* Using bacteria assembly accession names (note the --accession argument when using accession names instead of tax IDs):

```
zamp insilico -i zamp/data/bacteria-accs.txt \
-db greengenes2 --accession \
--fw-primer CCTACGGGNGGCWGCAG --rv-primer GACTACHVGGGTATCTAATCC 
```

* Using fungi tax IDs (requires additional ITS amplicon-specific parameters to adjust the amplicon size):

```
zamp insilico -i zamp/data/fungi-taxa.txt \
-db unite_db_v10 \ 
--fw-primer CYHRGYYATTTAGAGGWMSTAA --rv-primer RCKDYSTTCWTCRWYGHTGB \
--minlen 50 --maxlen 900
```

* Using a query term. In this example, 100 assemblies will be downloaded per taxon (`nb 100`) including non-reference assemblies (`not-only-ref`):

```
zamp insilico -i "lactobacillus" \
-db ezbiocloud \
--fw-primer CCTACGGGNGGCWGCAG --rv-primer GACTACHVGGGTATCTAATCC \
--replace-empty -nb 100 --not-only-ref
```


************************************************************************
Output
************************************************************************
The pipeline gathers information on available assemblies for the requested taxIDs in the `assembly_finder` folder.

The output of the in-silico amplification is in `Insilico` folder, and contains the following subfolders:

- PCR: contains the output of in-silico PCR amplification
- 2_denoised: output of clustering and denoising into representative sequences, and count tables
- 3_classified: output of taxonomic classification and tables comparing expected and obtained taxonomic assignations (`InSilico_compare_tax.tsv` and `InSilico_compare_tax_long.tsv`.)
