
.. _under_the_hood:


########################################################################
Under the hood
########################################################################

************************************************************************
Snakemake, environments and containers
************************************************************************



`Snakemake <https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html>`_ is the center-piece of this pipeline. 
Snakemake is a Python-based workflow-manager that enables the processing of a large set of amplicon-based metagenomics sequencing reads into actionable outputs. 
Each step is defined as a rule in which input/output files, software dependencies (Conda or containers), scripts and command-lines are specified (See `snakemake's docs <https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html>`_ for more details).




`Conda environments <https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html>`_

Conda is a language-independent package and environment management tool. A Conda environment is a collection of installed Conda packages. For example, a research project might require VSEARCH 2.20.0 and its dependencies, whereas another environment associated with a completed project might necessitate the use of VSEARCH 2.15. Changing the environment, has no effect on the others.
Switching between environments is simple because they can be easily activated or deactivated. 


`Apptainer containers <https://apptainer.org/docs/user/latest/>`_ 

The concept of reproducible analysis in bioinformatics extends beyond good documentation and code sharing. Analyses typically depend on an entire environment with numerous tools, libraries, and settings. 
Storage, reuse, and sharing environments via container software such as Docker and Singularity could improve reproducibility and productivity. 
By using containers (apptainer, docker, podman ...), users can create a single executable file that contains all aspects of their environment and allows to safely run environments from a variety of resources without requiring privileged access. 



************************************************************************
Logging and traceability
************************************************************************

logs
=======================================================================
Upon each execution, *zAMP* automatically creates a log file where all the standard output is recorded::

    zamp_out/zamp.log

config file
=======================================================================
In addition to logs, *zAMP* copies a config file listing all the parameters used in the run unde  ::

    zamp_out/config.yaml


************************************************************************
Sequencing reads QC
************************************************************************

QC rules assess the sequencing quality of all each sample with FastQC [1]_. Then, a MultiQC [2]_ report generates a report for each sequencing run (based on "run" column indicated  in *sample sheet* ). 
A global MultiQC report is generated as well, but without interactive features to deal with the high number of samples  


************************************************************************
Denoising
************************************************************************



Vsearch (OTU clustering)
=======================================================================

PANDAseq
-----------------------------------------------------------------------

Vsearch
-----------------------------------------------------------------------



DADA2 (ASV denoising)
=======================================================================

cutadapt
-----------------------------------------------------------------------

DADA2
-----------------------------------------------------------------------



************************************************************************
Taxonomic assignment
************************************************************************

reference database
=======================================================================

classifiers
=======================================================================



************************************************************************
Post-processing
************************************************************************


Taxonomic filtering
=======================================================================


Rarefaction
=======================================================================


Phylogenetic tree generation
=======================================================================


Taxonomic collapsing
=======================================================================


Normalization and abundance-based filtering
=======================================================================


Exports
=======================================================================


Fromatting
=======================================================================

Wide to long melting
-----------------------------------------------------------------------

transpose_and_meta_count_table
-----------------------------------------------------------------------

Qiime2 formats
-----------------------------------------------------------------------


************************************************************************
Picrust2
************************************************************************




************************************************************************
References
************************************************************************
.. [1] Andrews S, Krueger F, Seconds-Pichon A, Biggins F, Wingett S. FastQC. A quality control tool for high throughput sequence data. Babraham Bioinformatics. Babraham Institute. 2015. 
.. [2] Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: Summarize analysis results for multiple tools and samples in a single report. Bioinformatics. 2016; 




