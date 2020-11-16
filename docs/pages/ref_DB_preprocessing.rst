.. DB_preprocessing:
########################################################################
Taxonomic reference database preprocessing
########################################################################

************************************************************************
Rational:    
************************************************************************
After processing of the reads by a metagenomic pipeline, we expect amplicon sequences (OTUs or ASVs) to be assigned to the lowest possible taxonomic level (Species). However, all species cannot be differentiated without ambiguities based on short amplicon sequences. Thus, it is expected for some species to have more or less exactly the same sequence on the studied gene, and even more on the short fragment amplified and sequenced in amplicon-based metagenomics. Furthermore, the classifiers included in RSP4ABM (original RDP[1]_, RDP integrated in Qiime[2]_ and Decipher IDTAXA[3]_) all have specific formatting requirements or require a pre-training. Thus, to improve the taxonomic classification and to adapt the format of the provided reference database to multiple classifiers, RSP4ABM include a dedicated reference database preprocessing workflow. 

Working principle:
-----------------------------------------------------------------------

The user indicates to the pipeline the sequences of the used PCR primers and the path to the reference database fasta file and taxonomy annotation file to be formatted. 

With that and using tools from Qiime2[4]_ and VSEARCH[5]_, as well as home-made scripts, the pipeline will first extract the amplicon matching the used primes. Then, it will compare all the predicted. In cases where the exact same amplicon is predicted for multiple taxa, it unite their taxaonomy. For this, it will collapse together the genera/species 




Title3
=======================================================================

Title4
-----------------------------------------------------------------------

.. [1] Wang Q, Garrity GM, Tiedje JM, Cole JR. Naïve Bayesian classifier for rapid assignment of rRNA sequences into the new bacterial taxonomy. Appl Environ Microbiol. 2007. 
.. [2] Caporaso JG, Kuczynski J, Stombaugh J, Bittinger K, Bushman FD, Costello EK, et al. QIIME allows analysis of high-throughput community sequencing data. Nature Methods. 2010. 
.. [3] Murali A, Bhargava A, Wright ES. IDTAXA: A novel approach for accurate taxonomic classification of microbiome sequences. Microbiome. 2018.
.. [4] Estaki M, Jiang L, Bokulich NA, McDonald D, González A, Kosciolek T, et al. QIIME 2 Enables Comprehensive End-to-End Analysis of Diverse Microbiome Data and Comparative Studies with Publicly Available Data. Curr Protoc Bioinforma. 2020. 
.. [5] Rognes T, Flouri T, Nichols B, Quince C, Mahé F. VSEARCH: A versatile open source tool for metagenomics. PeerJ. 2016.