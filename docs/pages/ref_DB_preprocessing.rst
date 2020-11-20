
.. _DB_preprocessing:

########################################################################
Taxonomic reference database preprocessing
########################################################################

.. Note:: Provided command-line examples are given as examples and are valid for a standard unix bash terminal.


************************************************************************
Rational:    
************************************************************************
After processing of sequencing reads by a metagenomic pipeline, we expect amplicon sequences (OTUs or ASVs) to be assigned to the lowest possible taxonomic level (species). However, it is expected that different species might have a close to exact same sequence on the gene used as marker. Thus, all species cannot be differentiated without ambiguities based on marker genes, and that even more on the short fragment amplified and sequenced in amplicon-based metagenomics. 

The classifiers integrated in RSP4ABM (original *RDP* [1]_, *RDP* integrated in *QIIME* [2]_ and, *Decipher IDTAXA* [3]_) all have specific formatting requirements and the two last require an initial training. 

Thus, to improve the taxonomic classification and to provide to the end-user an information regarding the risk of confusion of certain taxa, a dedicated workflow of the pipeline will fuse taxa represented by identical sequences. This workflow also adapts the format of the database to make it compatible with multiple classifiers and conducts the training required by some of them. 

************************************************************************
Working principle:
************************************************************************

The user indicates to the pipeline: 

- the sequences of the used PCR primers and the path to the input reference database `fasta <https://en.wikipedia.org/wiki/FASTA_format>`_ file and taxonomy annotation file to be formatted.

Based on this information and using tools from *Cutadapt* [4]_ and *VSEARCH* [5]_, as well as home-made R [6]_ scripts, the pipeline first extracts the amplicon matching the used primes. Then, it unifies the taxonomy: in cases where the exact same amplicon is predicted for multiple taxa, it collapses together their identifiers at the genus/species (up to a user-defined number of occurrences) (#DAJ what happens when the threshold is reached?). An error is raised when the same sequence is observed across different families or ranks above (#DAJ what happens to those amplicons ? The pipline stops or they are discarded?).

In addition, the pipeline formats the database and executes the pre-training required for the original *RDP* classier as well as *Decipher IDTAXA*.

Finally, the  pipeline will make a copy of the original database as well as computes hashes of all files for traceability purposes. 

************************************************************************
Execution:
************************************************************************

A dedicated workflow is embedded in RSP4ABM for database preprocessing. This workflow must be run only one time for each set of PCR primer and reference database. 

First, the user must retrieve a database in the right format and a *config* file must be defined. Then, provided that the pipeline was properly setup (*see* :ref:`setup`), the dedicated workflow can be executed. 


Taxonomy database
=======================================================================

RSP4ABM requires a reference taxonomic database for classification. The database provided to RSP4ABM must be organized into two files, following the original *QIIME* format: 

- a `Reference sequences`_ `fasta <https://en.wikipedia.org/wiki/FASTA_format>`_ file.
- a `Reference taxonomy`_ text file describing with the taxonomic classification of these sequences.  


Reference sequences
-----------------------------------------------------------------------

The first file must be a `fasta file <https://en.wikipedia.org/wiki/FASTA_format>`_ with reference genomic sequences. The description of each sequence must be an unique sequence identifier.


*For instance*::

    >1
    CTGNCGGCGTGCCTAACACATNCAAGTCGAGCGGTGCTACGGAGGTCTTCGGACTGAAGTAGCATAGCGGCGGACGGGTGAGTAATACACAGGAACGTGCCCCTTGGAGGCGGATAGCTGTGGGAAACTGCAGGTAATCCGCCGTAAGCTCGGGAGAGGAAAGCCGGAAGGCGCCGAGGGAGCGGCCTGTGGCCCATCAGGTAGTTGGTAGGGTAAGAGCCTACCAAGCCGACGACGGGTAGCCGGTCTGAGAGGATGGACGGCCACAAGGGCACTGAGACACGGGCCCTACTCCTACGGGAGGCAGCAGTGGGGGATATTGGACAATGGGCGAAAGCCTGATCCAGCGACGCCGCGTGAGGGACGAAGTCCTTCGGGACGTAAACCTCTGTTGTAGGGGAAGAAGACAGTGACGGTACCCTACGAGGAAGCCCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGGNCGAGCGTTACCCGGAATCACTGGGCGTAAAGGGTGCGTA
    >2
    AACGAACGCTGGCGGCAGGCTTAACACATGCAAGTCGAACGAAGTCTTCGGACTTAGTGGCGCACGGGTGAGTAACACGTGGGAATATACCTCTTGGTGGGGAATAACGTCGGGAAACTGACGCTAATACCGCATACGCCCTTCGGGGGAAAGATTTATCGCCGAGAGATTAGCCCGCGTCCGATTAGCTAGTTGGTGAGGTAATGGCTCACCAAGGCGACGATCGGTAGCTGGTCTGAGAGGATGATCAGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGGACAATGGGCGAAAGCCTGATCCAGCCATGCCGCGTGAGTGATGAAGGCCTTAGGGTTGTAAAGCTCTTTCACCCACGACGATGATGACGGTAGTGGGAGAAGAAGCCCCGGCTAACTTCGTGCCAGCAGCCGCGGTAATACGAAGGGGGCTAGCGTTGTTCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTGGTCATAGTCAGAAGTGAAAGCCCTGGGCTCAACCCGGGAATTGCTTTTGATACTGGACCGCTAGAATCACGGAGAGGGTAGTGGAATTCCGAGTGTAGAGGTGAAATTCGTAGATATTCGGAAGAACACCAGTGGCGAAGG



Reference taxonomy
-----------------------------------------------------------------------

The second file must be a text file where the first column is the sequence identifier and the second represents its 7 levels of taxonomy separated by ";" (Kingdom;Phylum;Class;Order;Family;Genus;Genus Species). Both columns must be separated by a tabulation.

*For instance*::

    1   Bacteria;Proteobacteria;Alphaproteobacteria;Rhodospirillales;Rhodospirillaceae;Magnetospirillum;Magnetospirillum magnetotacticum
    2   Bacteria;Fusobacteria;Fusobacteria_c;Fusobacteriales;Fusobacteriaceae;Fusobacterium;Fusobacterium nucleatum



Find your own reference database
-----------------------------------------------------------------------
We do not provide a taxonomic reference database. However, here is a short, non-exhaustive, list of databases from which we could successfully prepare a database with our pipeline. 


*EzBioCloud (16S rRNA  - Bacteria)*

    `Website <https://www.ezbiocloud.net/resources/16s_download>`_ 

    `Publication <https://doi.org/10.1099/ijsem.0.001755>`_ 


*Silvia (16/18S rRNA, 23/28S rRNA - Bacteria and Eukarya )*

    `Website <https://www.arb-silva.de/download/arb-files/>`_ 

    `Publication <https://doi.org/10.1093/nar/gks1219>`_ 


*UNITE (ITS - Eukarya)* 

    `Website <https://unite.ut.ee/repository.php>`_ 


Working directory
=======================================================================

To execute the pipeline place yourself in any directory, but preferably not in the directory pipeline. It does not have to be where the input reference database files are, nor where you desire to save the output (these locations will be defined in the `config file`_ .) 

*for instance*::

    # create a directory to run the database preprocessing workflow
    $ mkdir DB_processing
    # Change into this new directory
    $ cd DB_processing


Config file
=======================================================================

As for the main pipeline (#DAJ not sure what is the main pipeline here? The genomic?), parameters must be provided in an *config file* in the *.yaml* format. Please adapt the following template to your situation.

*for instance*::

    # Open a graphic text editor and create a config file. Once opened, copy and paste the example below and adapt it. 
    $ gedit config_DB.yaml

    # Or use a command-line text editor, e.g. 
    # $ nano config_DB.yaml 



.. literalinclude:: ../../ressources/template_files/config_DB.yaml
    :language: yaml


Pipeline execution
=======================================================================

Once the reference database in the right format downloaded and the *config file* prepared, the database preprocessing pipeline can be executed. 


.. literalinclude:: ../../ressources/template_files/DB_snakemake_bash_command.sh 
    :language: bash



************************************************************************
Working with without pipeline preprocessing?:
************************************************************************

TO BE EXPLAINED!

************************************************************************
References:
************************************************************************

.. [1] Wang Q, Garrity GM, Tiedje JM, Cole JR. Naïve Bayesian classifier for rapid assignment of rRNA sequences into the new bacterial taxonomy. Appl Environ Microbiol. 2007. 
.. [2] Caporaso JG, Kuczynski J, Stombaugh J, Bittinger K, Bushman FD, Costello EK, et al. QIIME allows analysis of high-throughput community sequencing data. Nature Methods. 2010. 
.. [3] Murali A, Bhargava A, Wright ES. IDTAXA: A novel approach for accurate taxonomic classification of microbiome sequences. Microbiome. 2018.
.. [4] Compeau PEC, Pevzner PA, Tesler G, Papoutsoglou G, Roscito JG, Dahl A, et al. Cutadapt removes adapter sequences from high-throughput sequencing reads kenkyuhi hojokin gan rinsho kenkyu jigyo. EMBnet.journal. 2013.  
.. [5] Rognes T, Flouri T, Nichols B, Quince C, Mahé F. VSEARCH: A versatile open source tool for metagenomics. PeerJ. 2016
.. [6] R Core Development Team. R: A language and environment for statistical computing. Vienna, Austria. 2019. 
