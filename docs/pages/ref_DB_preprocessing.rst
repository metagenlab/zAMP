
.. _tax_DB:

########################################################################
Taxonomic reference databases
########################################################################


Database processing principle
---------------------------------

zAMP offers the possibility to train classifiers like *RDP* [1]_, *QIIME* [2]_ and *Decipher* [3]_ on full length sequences or subdomains like the V3-V4 region in the 16S rRNA gene.

Classification accuracy depends on classifiers parameters and the region on which they were trained on [4]_. In fact, Bokulich et al. [4]_ demonstrated that training classifiers on specific regions leads to enhanced accuracy compared to using the full length sequences (for short reads).

Therefore, by default, zAMP extracts primer amplified regions with cutadapt [5]_, dereplicates and clusters sequences with vseach [6]_, and adapts taxonomy according to these clusters.

Database processing flowchart:

.. mermaid::
    :align: center

    flowchart TD
    A["dna-sequences.fasta"]
    taxonomy["taxonomy.tsv"]
    taxonomy --> clean("Clean taxonomy")
    clean --> B
    A --> B{Extract region ?}
    B --> |Yes| C("Extract regions 
    (cutadapt)")
    B --> |No| D(Copy files)
    C --> E("Cluster
    (vsearch)")
    E --> F("Derep and merge taxonomy")
    F --> G("Train classifiers
    (RDP, Decipher...) ")

Custom database
------------

Taxonomy table in QIIME format::
    
        1   Bacteria;Proteobacteria;Alphaproteobacteria;Rhodospirillales;Rhodospirillaceae;Magnetospirillum;Magnetospirillum magnetotacticum
        2   Bacteria;Fusobacteria;Fusobacteria_c;Fusobacteriales;Fusobacteriaceae;Fusobacterium;Fusobacterium nucleatum
    
Fasta file::

        >1
        CTGNCGGCGTGCCTAACACATNCAAGTCGAGCGGTGCTACGGAGGTCTTCGGACTGAAGTAGCATAGCGGCGGACGGGTGAGTAATACACAGGAACGTGCCCCTTGGAGGCGGATAGCTGTGGGAAACTGCAGGTAATCCGCCGTAAGCTCGGGAGAGGAAAGCCGGAAGGCGCCGAGGGAGCGGCCTGTGGCCCATCAGGTAGTTGGTAGGGTAAGAGCCTACCAAGCCGACGACGGGTAGCCGGTCTGAGAGGATGGACGGCCACAAGGGCACTGAGACACGGGCCCTACTCCTACGGGAGGCAGCAGTGGGGGATATTGGACAATGGGCGAAAGCCTGATCCAGCGACGCCGCGTGAGGGACGAAGTCCTTCGGGACGTAAACCTCTGTTGTAGGGGAAGAAGACAGTGACGGTACCCTACGAGGAAGCCCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGGNCGAGCGTTACCCGGAATCACTGGGCGTAAAGGGTGCGTA
        >2
        AACGAACGCTGGCGGCAGGCTTAACACATGCAAGTCGAACGAAGTCTTCGGACTTAGTGGCGCACGGGTGAGTAACACGTGGGAATATACCTCTTGGTGGGGAATAACGTCGGGAAACTGACGCTAATACCGCATACGCCCTTCGGGGGAAAGATTTATCGCCGAGAGATTAGCCCGCGTCCGATTAGCTAGTTGGTGAGGTAATGGCTCACCAAGGCGACGATCGGTAGCTGGTCTGAGAGGATGATCAGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGGACAATGGGCGAAAGCCTGATCCAGCCATGCCGCGTGAGTGATGAAGGCCTTAGGGTTGTAAAGCTCTTTCACCCACGACGATGATGACGGTAGTGGGAGAAGAAGCCCCGGCTAACTTCGTGCCAGCAGCCGCGGTAATACGAAGGGGGCTAGCGTTGTTCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTGGTCATAGTCAGAAGTGAAAGCCCTGGGCTCAACCCGGGAATTGCTTTTGATACTGGACCGCTAGAATCACGGAGAGGGTAGTGGAATTCCGAGTGTAGAGGTGAAATTCGTAGATATTCGGAAGAACACCAGTGGCGAAGG

Default command::

    zamp db --taxonomy taxonomy.tsv \
            --fasta sequences.fasta \
            --fw-primer CCTACGGGNGGCWGCAG \
            --rv-primer GACTACHVGGGTATCTAATCC \
            -o processed-db

Skip primer amplified region extraction::

    zamp db --taxonomy taxonomy.tsv \
        --fasta sequences.fasta \
        --no-processing \
        -o unprocessed-db

Available databases
-------------------

.. Note:: Processed and unprocessed SILVA, Greengenes2 and UNITE will be made available soon

Here is a short, non-exhaustive, list of databases from which we could successfully prepare a database:

* EzBioCloud (16S rRNA  - Bacteria)

    `Website <https://www.ezbiocloud.net/resources/16s_download>`_ 

    `Publication <https://doi.org/10.1099/ijsem.0.001755>`_ 


* SILVA (16/18S rRNA, 23/28S rRNA - Bacteria and Eukarya )

    `Website <https://www.arb-silva.de/download/arb-files/>`_ 

    `Publication <https://doi.org/10.1093/nar/gks1219>`_ 


* UNITE (ITS - Eukarya)

    `Website <https://unite.ut.ee/repository.php>`_ 
    `Publication <https://doi.org/10.1093/nar/gkad1039>`_ 



* Eukaryome (ITS - Eukarya)

    `Website <https://eukaryome.org/download/>`_
    `Publication <https://doi.org/10.1093/database/baae043>`_ 

.. Note:: For Eukaryome, additional steps might be needed to filter a kingdom of interest (e.g. Fungi), or remove entries with incomplete taxonomy.

Parameters
---------------------------------

Command::

    zamp db -h

* The "--no-processing" parameter enables to skip the preprocessing and only format the provided database and train the classifiers. 
* "fw-primer" and "rv-primer" are fed to `cutadapt linked adapter argument <https://cutadapt.readthedocs.io/en/v3.0/guide.html#linked-adapters-combined-5-and-3-adapter>`_. 
* "--cutadapt_args_fw" and "--cutadapt_args_rv" allow to pass additional arguments to cutadapt, affecting the forward and reverse primer, respectively. It for instance allows to indicate which primer is optional <https://cutadapt.readthedocs.io/en/v3.0/guide.html#changing-which-adapters-are-required>`_. It is particularly useful when trying to extract ITS1 amplicons: the 5' universal primer is located on the SSU rRNA preceding the ITS region and thus is absent in ITS reference database. In this case, providing "--cutadapt_args_fw optional" enables to make it optional. 
* "errors" is fed to `cutadapt to define the number of accepted mismatches per primer <https://cutadapt.readthedocs.io/en/v3.0/guide.html#minimum-overlap-reducing-random-matches>`_. 
* "ampcov" is used with the length of the provided primers to feed `cutadapt with a minimal overlap <https://cutadapt.readthedocs.io/en/v3.0/guide.html#minimum-overlap-reducing-random-matches>`_.

Examples
---------

Bacteria
~~~~~~~~

**Greengenes2**

* Download::

    wget http://ftp.microbio.me/greengenes_release/2022.10/2022.10.backbone.full-length.fna.qza && \
    wget http://ftp.microbio.me/greengenes_release/2022.10/2022.10.backbone.tax.qza


* Decompress qza with `qiime2 export <https://docs.qiime2.org/2024.5/tutorials/exporting/)>`_::

    docker run -t -i -v $(pwd):/data quay.io/qiime2/tiny:2024.5 \
    qiime tools export \
    --input-path 2022.10.backbone.full-length.fna.qza \
    --output-path greengenes2 && \
    docker run -t -i -v $(pwd):/data quay.io/qiime2/tiny:2024.5 \
    qiime tools export \
    --input-path 2022.10.backbone.tax.qza --output-path greengenes2


* Prepare database::

    zamp db --fasta greengenes2/dna-sequences.fasta \
    --taxonomy greengenes2/taxonomy.tsv --name greengenes2 \
    --fw-primer CCTACGGGNGGCWGCAG --rv-primer GACTACHVGGGTATCTAATCC \
    -o greengenes2

Fungi
~~~~~~~~

**Unite ITS1**

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

**Eukaryome ITS2**

Similarly, to extract ITS2 from fungal databases such as Eukaryome, the reverse primer needs to be set as optional, because it is located on the LSU, which is absent in the database sequences::

    zamp db \
    --fasta QIIME2_EUK_ITS_v1.8.fasta \
    --taxonomy QIIME2_EUK_ITS_v1.8.txt \
    --name eukaryome \
    --fw-primer GCATCGATGAAGAACGCAGC --rv-primer TCCTCCGCTTATTGATATGC \
    --minlen 50 --maxlen 900 \
    --cutadapt_args_rv "optional" \
    -o eukaryome_ITS2


Output
------------------
Please, see <tax_DB_path>/<tax_DB_name>/QIIME/problematic_taxa.txt file for identical sequences that had taxonomic disagreeing identifiers above the genus rank. 


References
------------------

.. [1] Wang Q, Garrity GM, Tiedje JM, Cole JR. Naïve Bayesian classifier for rapid assignment of rRNA sequences into the new bacterial taxonomy. Appl Environ Microbiol. 2007. 
.. [2] Caporaso JG, Kuczynski J, Stombaugh J, Bittinger K, Bushman FD, Costello EK, et al. QIIME allows analysis of high-throughput community sequencing data. Nature Methods. 2010. 
.. [3] Murali A, Bhargava A, Wright ES. IDTAXA: A novel approach for accurate taxonomic classification of microbiome sequences. Microbiome. 2018.
.. [4] Bokulich, N. A., Kaehler, B. D., Rideout, J. R., Dillon, M., Bolyen, E., Knight, R., Huttley, G. A., & Gregory Caporaso, J. (2018). Optimizing taxonomic classification of marker-gene amplicon sequences with QIIME 2’s q2-feature-classifier plugin. In Microbiome (Vol. 6, Issue 1). Springer Science and Business Media LLC. https://doi.org/10.1186/s40168-018-0470-z 
.. [5] Compeau PEC, Pevzner PA, Tesler G, Papoutsoglou G, Roscito JG, Dahl A, et al. Cutadapt removes adapter sequences from high-throughput sequencing reads kenkyuhi hojokin gan rinsho kenkyu jigyo. EMBnet.journal. 2013.  
.. [6] Rognes T, Flouri T, Nichols B, Quince C, Mahé F. VSEARCH: A versatile open source tool for metagenomics. PeerJ. 2016

