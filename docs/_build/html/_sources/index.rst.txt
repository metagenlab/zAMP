.. RSP4ABM documentation master file, created by
   sphinx-quickstart on Fri May  8 17:11:11 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Reproducible Scalable Pipeline For Amplicon-based Metagenomics (RSP4ABM)
===============================================================================

RSP4ABM is a bioinformatic pipeline designed for convenient, reproducible and scalable amplicon-based metagenomics. RSP4ABM is a compilation of command-line tools and R packages that processes next-generation sequencing (NGS) reads into taxonomic profiles. Scripts launching these tools are embedded in a Snakemake pipeline for reproducible, scalable and easy "*one-command-line*" execution of all steps, from QC to basic plotting. Tools and packages are made available into Singularity containers and Conda environments for easy, reproducible management of dependencies. RST4ABM processes is capable of processing local *fastq* or web-stored sequences from the Sequence Reads Archive (SRA). Complementary workflows enables preprocessing of reference taxonomic databases as *in silico* validation of amplicon amplification and taxonomic assignment based on reference genomes. 

Please note that RSP4ABM compiles number of published tools. We would be very grateful that you acknowledge our work but also the contribution of the developers of these tools in your publications using RSP4ABM. 


.. toctree::
   :maxdepth: 2
   :glob:

   pages/under_the_hood.rst
   pages/setup.rst
   pages/ref_DB_preprocessing.rst
   pages/execution.rst
   pages/output.rst
   pages/downstream_analysis.rst
   pages/FAQ.rst

