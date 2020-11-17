.. DB_preprocessing:
########################################################################
Taxonomic reference database preprocessing
########################################################################

************************************************************************
Rational:    
************************************************************************
After processing of sequencing reads by a metagenomic pipeline, we expect amplicon sequences (OTUs or ASVs) to be assigned to the lowest possible taxonomic level (species). However, it is expected for some species to have more or less exactly the same sequence on the gene used as marker. Thus, all species cannot be differentiated without ambiguities based on marker genes, and that even more on the short fragment amplified and sequenced in amplicon-based metagenomics. 

The classifiers included in RSP4ABM (original RDP [1]_, RDP integrated in Qiime [2]_ and Decipher IDTAXA [3]_) all have specific formatting requirements and the two last require an initial training. 

Thus, to improve the taxonomic classification and to adapt the format of the provided reference database to these multiple classifiers, RSP4ABM includes a dedicated reference database preprocessing workflow. 

************************************************************************
Working principle:
************************************************************************

The user indicates to the pipeline the sequences of the used PCR primers and the path to the reference database fasta file and taxonomy annotation file to be formatted. 

With that and using tools from Qiime2 [4]_ and VSEARCH [5]_, as well as home-made scripts, the pipeline will first extract the amplicon matching the used primes. Then, it unifies the taxonomy: in cases where the exact same amplicon is predicted for multiple taxa, it collapses together their identifiers at the genus/species (up to a user-defined number of occurrences). An error is raised in cases where the same sequence is observed across different families. 

In addition, the pipeline formats the database and executes the pre-training required for the original RDP classier as well as Decipher IDTAXA.

Finally, the  pipeline will make a copy of the original database as well as computes hashes of all files for traceability purposes. 

************************************************************************
Execution:
************************************************************************

A dedicated workflow is embedded in RSP4ABM for database preprocessing. This workflow is to be run only one time for each set of PCR primer and reference database. It requires the pipeline to be properly :ref:`setup`  of a *config file* 

Config file
=======================================================================

.. literalinclude:: ../../ressources/template_files/config_DB.yaml
    :language: yaml





************************************************************************
References:
************************************************************************

.. [1] Wang Q, Garrity GM, Tiedje JM, Cole JR. Naïve Bayesian classifier for rapid assignment of rRNA sequences into the new bacterial taxonomy. Appl Environ Microbiol. 2007. 
.. [2] Caporaso JG, Kuczynski J, Stombaugh J, Bittinger K, Bushman FD, Costello EK, et al. QIIME allows analysis of high-throughput community sequencing data. Nature Methods. 2010. 
.. [3] Murali A, Bhargava A, Wright ES. IDTAXA: A novel approach for accurate taxonomic classification of microbiome sequences. Microbiome. 2018.
.. [4] Estaki M, Jiang L, Bokulich NA, McDonald D, González A, Kosciolek T, et al. QIIME 2 Enables Comprehensive End-to-End Analysis of Diverse Microbiome Data and Comparative Studies with Publicly Available Data. Curr Protoc Bioinforma. 2020. 
.. [5] Rognes T, Flouri T, Nichols B, Quince C, Mahé F. VSEARCH: A versatile open source tool for metagenomics. PeerJ. 2016