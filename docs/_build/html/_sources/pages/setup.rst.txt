########################################################################
Setup and system requirements
########################################################################

.. Note:: Provided command-line examples are given as example for a standard unix terminal in bash.

************************************************************************
Operating system and system resource 
************************************************************************

Operating system
=======================================================================
RST4ABM was designed on *Ubuntu 18.04** (Linux) but should compatible with all system capable of installing the dependencies listed on this page.

RAM memory
=======================================================================
Some tools embedded in RST4ABM can be quite demanding regarding RAM memory. The actual requirement depends on your dataset and is influenced parameters indicated in the config file. The bottleneck usually is the taxonomic assignment. Factors which can raise your RAM requirements are:
- the number of samples
- the bacterial diversity within your samples
- the number of cores (start by reducing the number of core in snakemake if facing issues)

.. Hint:: In practice, 16 to 32 GB are usually required. 


************************************************************************
Software dependencies
************************************************************************

*Git*
=======================================================================

What for?
-----------------------------------------------------------------------

*Git* is required to download (clone) RSP4ABM. 


Install
-----------------------------------------------------------------------

*Git* is available by default in operating systems. 

*Test if it is installed*::

    #To test if git is installed, make it print its version. will fail if not installed
    git --version

If it is not installed, `follow the indications on the git installation page. <https://git-scm.com/downloads>`_



*Conda*
=======================================================================

What for?
-----------------------------------------------------------------------

`Conda is a convenient python-based package and environment manager. <https://docs.conda.io/en/latest>`_
It  enables the easy installation of *Snakemake*. Furthermore, it can be used (as an alternative to *Singularity containers*) by *Snakemake* and RSP4ABM to install all the packages required to the execution of the pipeline.


Install
-----------------------------------------------------------------------
`Follow the Miniconda3 installation recommendation <https://docs.conda.io/en/latest/miniconda.html>`_



*Snakemake*
=======================================================================

What for?
-----------------------------------------------------------------------
*RSP4ABM* is a *Snakemake* pipeline. Then, is must be installed an available for each execution of the pipeline. 


Install
-----------------------------------------------------------------------
Follow indications on Snakemake `installation page <https://snakemake.readthedocs.io/en/stable/getting_started/installation.html>`_




*Reference database*

RSP4ABM requires a reference taxonomic database for taxonomic classification. The original reference database provided to RSP4ABM must be organized into two files: 

Reference sequences
-----------------------------------------------------------------------

The first file must be a `*.fasta* file <https://en.wikipedia.org/wiki/FASTA_format>` with reference genomic sequences. The description of each sequence must be an unique sequence identifier.


*For instance*::

    >1
    CTGNCGGCGTGCCTAACACATNCAAGTCGAGCGGTGCTACGGAGGTCTTCGGACTGAAGTAGCATAGCGGCGGACGGGTGAGTAATACACAGGAACGTGCCCCTTGGAGGCGGATAGCTGTGGGAAACTGCAGGTAATCCGCCGTAAGCTCGGGAGAGGAAAGCCGGAAGGCGCCGAGGGAGCGGCCTGTGGCCCATCAGGTAGTTGGTAGGGTAAGAGCCTACCAAGCCGACGACGGGTAGCCGGTCTGAGAGGATGGACGGCCACAAGGGCACTGAGACACGGGCCCTACTCCTACGGGAGGCAGCAGTGGGGGATATTGGACAATGGGCGAAAGCCTGATCCAGCGACGCCGCGTGAGGGACGAAGTCCTTCGGGACGTAAACCTCTGTTGTAGGGGAAGAAGACAGTGACGGTACCCTACGAGGAAGCCCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGGNCGAGCGTTACCCGGAATCACTGGGCGTAAAGGGTGCGTA
    >2
    AACGAACGCTGGCGGCAGGCTTAACACATGCAAGTCGAACGAAGTCTTCGGACTTAGTGGCGCACGGGTGAGTAACACGTGGGAATATACCTCTTGGTGGGGAATAACGTCGGGAAACTGACGCTAATACCGCATACGCCCTTCGGGGGAAAGATTTATCGCCGAGAGATTAGCCCGCGTCCGATTAGCTAGTTGGTGAGGTAATGGCTCACCAAGGCGACGATCGGTAGCTGGTCTGAGAGGATGATCAGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGGACAATGGGCGAAAGCCTGATCCAGCCATGCCGCGTGAGTGATGAAGGCCTTAGGGTTGTAAAGCTCTTTCACCCACGACGATGATGACGGTAGTGGGAGAAGAAGCCCCGGCTAACTTCGTGCCAGCAGCCGCGGTAATACGAAGGGGGCTAGCGTTGTTCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTGGTCATAGTCAGAAGTGAAAGCCCTGGGCTCAACCCGGGAATTGCTTTTGATACTGGACCGCTAGAATCACGGAGAGGGTAGTGGAATTCCGAGTGTAGAGGTGAAATTCGTAGATATTCGGAAGAACACCAGTGGCGAAGG

Reference taxonomy
-----------------------------------------------------------------------

The second file must be a text file where the first column is the sequence identifier and the second must be 7 levels of taxonomy separated by ";" (Kingdom;Phylum;Class;Order;Family;Genus;Genus Species). Both columns must be separated by a tabulation.

*For instance*::

    1   Bacteria;Proteobacteria;Alphaproteobacteria;Rhodospirillales;Rhodospirillaceae;Magnetospirillum;Magnetospirillum magnetotacticum
    2   Bacteria;Fusobacteria;Fusobacteria_c;Fusobacteriales;Fusobacteriaceae;Fusobacterium;Fusobacterium nucleatu