
.. _DB_preprocessing:

########################################################################
Taxonomic reference database preprocessing
########################################################################

.. Note:: Provided command-line examples are given as examples and are valid for a standard unix bash terminal.

.. Hint:: In principle, the preprocessing of the reference database exposed here will have to be conducted only once. 

************************************************************************
Rational:    
************************************************************************
After processing of sequencing reads by a metagenomic pipeline, we expect amplicon sequences (OTUs or ASVs) to be assigned to the lowest possible taxonomic level (species). However, it is expected that different species might have close to the exact same sequence on the gene used as marker. Thus, all species cannot be differentiated without ambiguities based on marker genes, and that even more on the short fragment amplified and sequenced in amplicon-based metagenomics. 

The classifiers integrated in RSP4ABM (original *RDP* [1]_, *RDP* integrated in *QIIME* [2]_ and, *Decipher IDTAXA* [3]_) all have specific formatting requirements and the two last require an initial training. 

Thus, to improve the taxonomic classification and to provide to the end-user an information regarding the risk of confusion of certain taxa, a dedicated workflow of the pipeline will fuse taxa represented by identical sequences. This workflow also adapts the format of the database to make it compatible with multiple classifiers and conducts the training required by some of them. 

************************************************************************
Working principle:
************************************************************************

The user indicates to the pipeline: 

- the sequences of the used PCR primers and the path to the input reference database `fasta <https://en.wikipedia.org/wiki/FASTA_format>`_ file and taxonomy annotation file to be formatted.

Based on this information and using tools from *Cutadapt* [4]_ and *VSEARCH* [5]_, as well as home-made R [6]_ scripts, the pipeline first extracts the amplicon matching the used primes. Then, it unifies the taxonomy: in cases where the exact same amplicon is predicted for multiple taxa, it collapses together their identifiers at the genus/species. Above a certain number of genus/species (defined by the user when preprocessing de database), the taxonomic identifier is remplaced by a `placeholder <https://en.wikipedia.org/wiki/Placeholder_name>`_ ("gen."/"sp."). In all cases, the number of taxonomic identifers is indicated as well between parentheses. 


Cases where identical sequences belong to different families or above after collapsing of identifiers are writted in a dedicated file. Futhermore, all ranks below the rank of disagreement are remplaced by an indicator of the disagreement. For instance::

    # A exemple of a sequence found in two distinct families (Sphingomonadaceae/Erythrobacteraceae). "Disc.Fam_Sphingomonadaceae/Erythrobacteraceae(2)" at the genus and species levels indicate this discrepancy.

    Bacteria;Proteobacteria;Alphaproteobacteria;Sphingomonadales;Sphingomonadaceae/Erythrobacteraceae;Disc.Fam_Sphingomonadaceae/Erythrobacteraceae(2);Disc.Fam_Sphingomonadaceae/Erythrobacteraceae(2)

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

To execute the pipeline place yourself in any directory. It does not have to be where the input reference database files are, nor where you desire to save the output (these locations will be defined in the `config file`_.) 

*for instance*::

    # create a directory to run the database preprocessing workflow
    $ mkdir DB_processing
    # Change into this new directory
    $ cd DB_processing

.. Hint:: For traceability and reproducibility, create this working directory and place your processed taxonomy database in a location where it will not be erased by error.

Config file
=======================================================================

As for the :ref:`main workflow <pipeline_execution>`, parameters must be provided in an *config file* in the *.yaml* format. Please adapt the following template to your situation.


*for instance*::

    # Open a graphic text editor and create a config file. Once opened, copy and paste the example below and adapt it. 
    $ gedit config_DB.yaml

    # Or use a command-line text editor, e.g. 
    # $ nano config_DB.yaml 



.. literalinclude:: ../../ressources/template_files/config_DB.yaml
    :language: yaml


.. Hint:: The "extract_and_merge" parameter in config enable to skip the preprocessing and only to format the provided database and train the classifiers. 

.. Hint:: "forward_primer" and "reverse_primer" are fed to `cutadapt linked adapter argument <https://cutadapt.readthedocs.io/en/v3.0/guide.html#linked-adapters-combined-5-and-3-adapter>`_. It for instance allows to indicate which primer is optional <https://cutadapt.readthedocs.io/en/v3.0/guide.html#changing-which-adapters-are-required>`_. It is particularly useful when trying to extract V1V2 amplicons: the 5' primer can be located before the 16S rRNA sequence provided in reference database. In this case, providing "FPRIMERSEQUENCE;optional" to the "forward_primer" enables to make it optional. 

.. Hint:: "excepted_errors" is fed to `cutadapt to define the number of accepted mismatches per primer <https://cutadapt.readthedocs.io/en/v3.0/guide.html#minimum-overlap-reducing-random-matches>`_. The "amplicon_min_coverage" is used with the length of the provided primers to feed `cutadapt with a minimal overlap <https://cutadapt.readthedocs.io/en/v3.0/guide.html#minimum-overlap-reducing-random-matches>`_. Emprically, 4 for the "excepted_errors" and 0.9 for the "amplicon_min_coverage" seems to be good values, not to loose to much sequences. Yet that may be something to play with. 



Pipeline execution
=======================================================================

Once the reference database in the right format downloaded and the *config file* prepared, the database preprocessing pipeline can be executed. 


.. literalinclude:: ../../ressources/template_files/DB_snakemake_bash_command.sh 
    :language: bash



Observe output
=======================================================================
Based on what was indicated in the *config file*, the preprocessed database will be located in::

    <tax_DB_path>/<tax_DB_name>/

.. attention:: Please, observe the <tax_DB_path>/<tax_DB_name>/QIIME/problematic_taxa.txt file for identical sequences that had taxonomic disagreeing identifiers above the genus rank. 


Generated output
=======================================================================

::

    ├── dada2rdp # DB formatted for RDP implemented in DADA2
    │   ├── DADA2_DB_amp_taxonomy_Genus_species.txt
    │   ├── DADA2_DB_amp_taxonomy_King_to_Genus.txt
    │   ├── DADA2_DB_amp_taxonomy_King_to_Species.txt
    │   └── DADA2_DB.hash
    ├── decipher # DB formatted for IDTAXA (decipher)
    │   └── Decipher_DB_amp_taxonomy.fasta
    ├── logs # Logs of the DB preprocessing 
    │   ├── 2020
    │   │   └── 11
    │   │       └── 19
    │   │           └── 14_15_15_98
    │   │               ├── cmd.txt
    │   │               ├── config.yaml       
    │   │               ├── git.txt
    │   │               └── user.txt
    │   ├── dada2rdp
    │   │   └── DB_amp_taxonomy_dada2_prep.log
    │   ├── decipher
    │   │   ├── DB_amp_taxonomy_decipher_fasta.log
    │   │   └── DB_amp_taxonomy_decipher_tax_tree.log
    │   ├── QIIME
    │   │   ├── DB_cutadapt.txt
    │   │   ├── derep_and_merge.log
    │   │   └── vsearch_dereplicate_ampli.log
    │   └── RDP
    │       ├── formatted_tax_table.log
    │       └── RDP_train.log
    ├── master # Copies of original DB for traceability
    │   ├── original.hash
    │   ├── original_seqs.fasta
    │   └── original_tax.txt
    ├── QIIME # Processed DB in original QIIME format
    │   ├── DB_amp_all_taxonomy.txt # All identifiers collapsed, without placeholders
    │   ├── DB_amp.fasta # deplicated sequences after amplicon extraction
    │   ├── DB_amp_taxonomy.txt # final taxonomy after collapsing 
    │   ├── DB_amp.uc # vsearch based clustering of identical sequences
    │   ├── DB_formatted.hash 
    │   ├── dna-sequences.fasta # extracted amplicons
    │   └── problematic_taxa.txt # sequences with collapsed taxa above the genus rank
    └── RDP # DB formatted for the original RDP
        ├── bergeyTrainingTree.xml
        ├── formatted_tax_table.tsv
        ├── genus_wordConditionalProbList.txt
        ├── logWordPrior.txt
        ├── RDP_DB.hash
        ├── ready4train_lineages.txt
        ├── ready4train_seqs.fasta
        ├── rRNAClassifier.properties
        └── wordConditionalProbIndexArr.txt



************************************************************************
References:
************************************************************************

.. [1] Wang Q, Garrity GM, Tiedje JM, Cole JR. Naïve Bayesian classifier for rapid assignment of rRNA sequences into the new bacterial taxonomy. Appl Environ Microbiol. 2007. 
.. [2] Caporaso JG, Kuczynski J, Stombaugh J, Bittinger K, Bushman FD, Costello EK, et al. QIIME allows analysis of high-throughput community sequencing data. Nature Methods. 2010. 
.. [3] Murali A, Bhargava A, Wright ES. IDTAXA: A novel approach for accurate taxonomic classification of microbiome sequences. Microbiome. 2018.
.. [4] Compeau PEC, Pevzner PA, Tesler G, Papoutsoglou G, Roscito JG, Dahl A, et al. Cutadapt removes adapter sequences from high-throughput sequencing reads kenkyuhi hojokin gan rinsho kenkyu jigyo. EMBnet.journal. 2013.  
.. [5] Rognes T, Flouri T, Nichols B, Quince C, Mahé F. VSEARCH: A versatile open source tool for metagenomics. PeerJ. 2016
.. [6] R Core Development Team. R: A language and environment for statistical computing. Vienna, Austria. 2019. 
