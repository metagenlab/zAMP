
.. _DB_preprocessing:

########################################################################
Taxonomic reference database preprocessing
########################################################################

.. Note:: Provided command-line examples are given as examples and are valid for a standard unix bash terminal.

.. Hint:: In principle, the preprocessing of the reference database exposed here will have to be conducted only once. 

************************************************************************
Rational:    
************************************************************************
After processing sequencing reads by a metagenomic pipeline, we expect amplicon sequences (OTUs or ASVs) to be assigned to the lowest possible taxonomic level (species). However, it is expected that different species might have close to the exact same sequence on the gene used as a marker. Thus, all species cannot be differentiated without ambiguities based on marker genes, and that is even more on the short fragment amplified and sequenced in amplicon-based metagenomics. 

The classifiers integrated into zAMP (original *RDP* [1]_, *RDP* integrated into *QIIME* [2]_ and, *Decipher IDTAXA* [3]_) all have specific formatting requirements and the two last require initial training. 

Thus, to improve the taxonomic classification and to provide the end-user an information regarding the risk of confusion of certain taxa, a dedicated workflow of the pipeline will fuse taxa represented by identical sequences. This workflow also adapts the format of the database to make it compatible with multiple classifiers and conducts the training required by some of them. 

************************************************************************
Working principle:
************************************************************************

The user indicates to the pipeline: 

- the sequences of the used PCR primers and the path to the input reference database `fasta <https://en.wikipedia.org/wiki/FASTA_format>`_ file and taxonomy annotation file to be formatted.

Based on this information and using tools from *Cutadapt* [4]_ and *VSEARCH* [5]_, as well as homemade R [6]_ scripts, the pipeline first extracts the amplicon matching the used primes. Then, it unifies the taxonomy: in cases where the exact same amplicon is predicted for multiple taxa, it collapses together their identifiers at the genus/species. Above a certain number of genus/species (defined by the user when preprocessing the database), the taxonomic identifier will be replaced by a `placeholder <https://en.wikipedia.org/wiki/Placeholder_name>`_ ("gen."/"sp."). In all cases, the number of taxonomic identifiers is indicated as well between parentheses. 


Cases, where identical sequences belong to different families or above after collapsing of identifiers, are written in a dedicated file. Furthermore, all ranks below the rank of disagreement are replaced by an indicator of the disagreement. For instance::

    # An example of a sequence found in two distinct families (Sphingomonadaceae/Erythrobacteraceae). "Disc.Fam_Sphingomonadaceae/Erythrobacteraceae(2)" at the genus and species levels indicate this discrepancy.

    Bacteria;Proteobacteria;Alphaproteobacteria;Sphingomonadales;Sphingomonadaceae/Erythrobacteraceae;Disc.Fam_Sphingomonadaceae/Erythrobacteraceae(2);Disc.Fam_Sphingomonadaceae/Erythrobacteraceae(2)

In addition, the pipeline formats the database and executes the pre-training required for the original *RDP* classier as well as *Decipher IDTAXA*.

Finally, the  pipeline will make a copy of the original database as well as compute hashes of all files for traceability purposes. 

************************************************************************
Execution:
************************************************************************

A dedicated workflow is embedded in zAMP for database preprocessing. This workflow must be run only once for each set of PCR primers and reference database. 

First, the user must retrieve a database in the right format.


Taxonomy database
=======================================================================

zAMP requires a reference taxonomic database for classification. The database provided to zAMP must be organized into two files, following the original *QIIME* format: 

- a `Reference sequences`_ `fasta <https://en.wikipedia.org/wiki/FASTA_format>`_ file.
- a `Reference taxonomy`_ text file describing the taxonomic classification of these sequences.  


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

Unite proposes a release conveniently in QIIME format, ready to use for `zAMP db`.

*Eukaryome (ITS - Eukarya)*

    `Website <https://eukaryome.org/download/>`_

Eukaryome proposes a release conveniently in QIIME format, ready to use for `zAMP db`.


Pipeline execution
=======================================================================

The module is executed as `zamp db` and has the following arguments::
    Usage: zamp db [OPTIONS] [SNAKE_ARGS]...

        Prepare database files for zAMP

        Options:
        --fasta PATH                    Path to database fasta file  [required]
        --taxonomy PATH                 Path to tab seperated taxonomy file in QIIME
                                        format  [required]
        --name TEXT                     Comma seperated list of database names
                                        [required]
        --processing / --no-processing  Extract amplicon regions and merge taxonomy
                                        [default: processing]
        --tax-collapse TEXT             Dictionary of number of ranks to print limit
                                        when collapsing names   [default:
                                        {"Species": 5, "Genus": 6}]
        --fw-primer TEXT                Forward primer sequence to extract amplicon
                                        [required]
        --rv-primer TEXT                Reverse primer sequence to extract amplicon
                                        [required]
        --minlen INTEGER                Minimum amplicon length  [default: 300]
        --maxlen INTEGER                Maximum amplicon length  [default: 500]
        --ampcov FLOAT                  Minimum amplicon coverage  [default: 0.9]
        --errors FLOAT                  Maximum number of accepted primer
                                        mismatches, or float between 0 and 1
                                        [default: 0.1]
        --cutadapt_args_fw TEXT         Additional cutadapt arguments for forward
                                        primer
        --cutadapt_args_rv TEXT         Additional cutadapt arguments for reverse
                                        primer
        --rdp-mem TEXT                  Maximum RAM for RDP training  [default: 30g]
        --classifier [rdp|qiimerdp|dada2rdp|decipher]
                                        Which classifiers to train on the database
                                        [default: rdp, qiimerdp, dada2rdp]
        -o, --output PATH               Output directory  [default: zamp_out]
        --configfile TEXT               Custom config file [default:
                                        (outputDir)/config.yaml]
        -t, --threads INTEGER           Number of threads to use  [default: 1]
        --use-singularity / --no-use-singularity
                                        Use singularity containers for Snakemake
                                        rules  [default: use-singularity]
        --singularity-prefix PATH       Custom singularity container directory
        --use-conda / --no-use-conda    Use conda for Snakemake rules  [default: no-
                                        use-conda]
        --conda-prefix PATH             Custom conda env directory
        --snake-default TEXT            Customise Snakemake runtime args
        -h, --help                      Show this message and exit.


* The "processing / --no-processing" parameter in config enables to skip the preprocessing and only format the provided database and train the classifiers. 
* "fw-primer" and "rv-primer" are fed to `cutadapt linked adapter argument <https://cutadapt.readthedocs.io/en/v3.0/guide.html#linked-adapters-combined-5-and-3-adapter>`_. 
* "--cutadapt_args_fw" and "--cutadapt_args_rv" allow to pass additional arguments to cutadapt, affecting the forward and reverse primer, respectively. It for instance allows to indicate which primer is optional <https://cutadapt.readthedocs.io/en/v3.0/guide.html#changing-which-adapters-are-required>`_. It is particularly useful when trying to extract ITS1 amplicons: the 5' universal primer is located on the SSU rRNA preceding the ITS region and thus is absent in ITS reference database. In this case, providing "--cutadapt_args_fw optional" enables to make it optional. 
* "errors" is fed to `cutadapt to define the number of accepted mismatches per primer <https://cutadapt.readthedocs.io/en/v3.0/guide.html#minimum-overlap-reducing-random-matches>`_. 
* "ampcov" is used with the length of the provided primers to feed `cutadapt with a minimal overlap <https://cutadapt.readthedocs.io/en/v3.0/guide.html#minimum-overlap-reducing-random-matches>`_.


Usage cases
=======================================================================

bacteria usage cases...

**Unite ITS1 database**

Fungal ITS databases Unite v10 and Eukaryome v1.8 do not contain the adjacent SSU/LSU sequences (they contain 5.8S), where some of the commonly used PCR primers lie on. 
It is important to adjust the cutadapt parameters so that only the absent primer is optional.
In the following example, we prepare a database for fungal ITS1 from Unite Db. 
In this case, the forward primer (lying of the 18S) will not be present in most sequences of Unite/Eukaryome (but the reverse primer lying on the 5.8S is present); therefore we set the forward primer as optional; the extracted sequences will start at the available 5' of the database and end at the reverse primer::

    zamp db \
    --fasta sh_refs_qiime_unite_ver10_dynamic_04.04.2024.fasta \
    --taxonomy sh_taxonomy_qiime_unite_ver10_dynamic_04.04.2024.txt \
    --name unite \
    --fw-primer CYHRGYYATTTAGAGGWMSTAA --rv-primer RCKDYSTTCWTCRWYGHTGB \
    --minlen 50 --maxlen 900 \
    --cutadapt_args_fw "optional" \
    -o unite_ITS1

**Eukaryome ITS2 database**

Similarly, to extract ITS2 from fungal databases such as Eukaryome, the reverse primer needs to be set as optional, because it is located on the LSU, which is absent in the database sequences::
    zamp db \
    --fasta QIIME2_EUK_ITS_v1.8.fasta \
    --taxonomy QIIME2_EUK_ITS_v1.8.txt \
    --name eukaryome \
    --fw-primer GCATCGATGAAGAACGCAGC --rv-primer TCCTCCGCTTATTGATATGC \
    --minlen 50 --maxlen 900 \
    --cutadapt_args_rv "optional" \
    -o eukaryome_ITS2



Generated output
=======================================================================
Please, observe the <tax_DB_path>/<tax_DB_name>/QIIME/problematic_taxa.txt file for identical sequences that had taxonomic disagreeing identifiers above the genus rank. 


************************************************************************
References:
************************************************************************

.. [1] Wang Q, Garrity GM, Tiedje JM, Cole JR. Naïve Bayesian classifier for rapid assignment of rRNA sequences into the new bacterial taxonomy. Appl Environ Microbiol. 2007. 
.. [2] Caporaso JG, Kuczynski J, Stombaugh J, Bittinger K, Bushman FD, Costello EK, et al. QIIME allows analysis of high-throughput community sequencing data. Nature Methods. 2010. 
.. [3] Murali A, Bhargava A, Wright ES. IDTAXA: A novel approach for accurate taxonomic classification of microbiome sequences. Microbiome. 2018.
.. [4] Compeau PEC, Pevzner PA, Tesler G, Papoutsoglou G, Roscito JG, Dahl A, et al. Cutadapt removes adapter sequences from high-throughput sequencing reads kenkyuhi hojokin gan rinsho kenkyu jigyo. EMBnet.journal. 2013.  
.. [5] Rognes T, Flouri T, Nichols B, Quince C, Mahé F. VSEARCH: A versatile open source tool for metagenomics. PeerJ. 2016
.. [6] R Core Development Team. R: A language and environment for statistical computing. Vienna, Austria. 2019. 
