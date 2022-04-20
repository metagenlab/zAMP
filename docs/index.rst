.. RSP4ABM documentation master file, created by
   sphinx-quickstart on Fri May  8 17:11:11 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Reproducible Scalable Pipeline For Amplicon-based Metagenomics (RSP4ABM)
===============================================================================

RSP4ABM is a bioinformatic pipeline designed for convenient, reproducible and scalable amplicon-based metagenomics. RSP4ABM is a compilation of command-line tools and R packages that processes next-generation sequencing (NGS) from amplicon-based metagenomics reads into taxonomic profiles. Scripts running these tools are embedded in a `Snakemake <https://snakemake.readthedocs.io/en/stable/>`_ pipeline for reproducible, scalable and easy "*one-command-line*" execution of all steps, from QC to basic plotting. Tools and packages are made available into `Singularity <https://sylabs.io/guides/latest/user-guide/>`_ containers and `Conda <https://docs.conda.io/en/latest/>`_ environments for easy and reproducible management of dependencies. RST4ABM is capable of processing local sequencing reads as well as reads from the `Sequence Reads Archive (SRA) <https://www.ncbi.nlm.nih.gov/sra>`_. 

A complementary workflows enables preprocessing of reference taxonomic databases needed to classify taxa to amplicon sequences. Another workflow allows to predict *in silico* the expected taxonomic classification from selected reference genomes.

Please note that RSP4ABM is a compilation of a number of published tools. We would be very grateful that you acknowledge our work but also the contribution of the developers of these tools in your publications using RSP4ABM. 


.. toctree::
   :maxdepth: 2
   :glob:

   pages/setup.rst
   pages/ref_DB_preprocessing.rst
   pages/execution.rst
   pages/under_the_hood.rst
   pages/output.rst
   pages/downstream_analysis.rst
   pages/FAQ.rst
   pages/insilico_validation.rst

