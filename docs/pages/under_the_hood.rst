
.. _under_the_hood:
########################################################################
Under the hood
########################################################################

************************************************************************
Snakemake, environments and containers
************************************************************************

`**Snakemake**`<https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html<

Snakemake is the center-piece of this pipeline. Snakemakie is a Python-based worflow-manager that enables processing of a large set of amplicon-based metagenomics sequencing reads into actionable output. It relies on a system of rules which are as many required processing steps. Each rule specifies input files, a Conda environment (or as an alternative a Singularity container) that includes all required softwares, and command line or a script to be executed the expected output files.



:github:`Conda environments <envs>` 

Conda is a language-independent package and environment management tool. The Conda environment is a collection of installed Conda packages. For example, a research project might require VSEARCH 2.20.0 and its dependencies, whereas another environment associated with a completed project might necessitate the use of VSEARCH 2.15. Changing the environment, has no effect on the others.
Switching between environments is simple because they can be easily activated or deactivated. 


:github:`Singularity containers <envs/singularity>` 

The concept of reproducible analysis in bioinformatics extends beyond good documentation and code sharing. Analyses typically depend on an entire environment with numerous tools, libraries, and settings. Storage, reuse, and sharing environments via container software such as Docker and Singularity could improve reproducibility and productivity. By using singularity, users can create a single executable file that contains all aspects of their environment and allows to safely run environments from a variety of resources without requiring privileged access. 


************************************************************************
Input output definition 
************************************************************************

Input
=======================================================================

:github:`A Python script <rules/0_preprocessing/scripts/make_input_lists.py>` is run at each execution of the RSP4ABM to match sequencing read files to each sample. 

As exposed in :ref:`the sample selection section of the Pipeline execution page <sample_selection>`, this function can accept either local or SRA-hosted read files. 

Local reads
-----------------------------------------------------------------------

Reads are considered to be located locally :github:`by the script pointing toward the input read files <rules/0_preprocessing/scripts/make_input_list.py>` in presence of the "local_samples:" argument in *the config file*. This parameters must point to a :github:`sample sheet <ressources/template_files/example_local_samples.tsv>`, which is a spreadsheet in `tabulation-separated values (tsv) format <https://en.wikipedia.org/wiki/Tab-separated_values>`_ listing all the files in the analysis. The script matches by default the values found in the leftmost "Sample" column of the *sample sheet* with the filenames of the ".fastq.gz" locales in a directory defined by the "link_directory" parameter of the *config file*. This default behavior can be altered by two conditions encoded in the :github:`Python script <rules/0_preprocessing/scripts/make_input_list.py>` :

- In presence of a column named "OldSampleName", the match with the filenames of the sequencing read files is done with this column, instead of the leftmost "Sample" column. 
- In presence of a "R1" column, the absolute path to the forward reads in considered instead of a name-based matching of the sequencing read files. In this case, and for paired-end reads, a "R2" column must point to the reverse read files. 

In all cases (i.e. *Sample column match*, *OldSampleName column match* or *absolute paths indicated by R1 and R2 columns*), :github:`Snakemake rules <rules/0_preprocessing/get_inputs.rules>` will temporarily copy the read files into a *raw_reads* directory. 

SRA reads
-----------------------------------------------------------------------

Reads are considered by :github:`the Python scripts <rules/0_preprocessing/scripts/make_input_lists.py>` to be located on the Sequence Read Archive (`SRA <https://en.wikipedia.org/wiki/Sequence_Read_Archive>`_) in presence of the "sra_samples:" argument in the *config file*.

In this case :github:`Snakemake rules <rules/0_preprocessing/get_inputs.rules>` will use `SRA Toolkit <https://github.com/ncbi/sra-tools>`_ to download the reads and convert them to ".fastq.gz" format into the *raw_reads* directory.

Output
=======================================================================

Upon each execution of the pipeline, :github:`Python scripts <rules/0_preprocessing/scripts/make_output_lists.py>` will parse the *config file* and the *sample sheet* to generate lists of outputs. These lists are then fed to the :github:`pipeline Snakefile <Snakefile>` which instructs the pipeline the output to generate.


************************************************************************
Logging and traceability
************************************************************************

Snakemake logs
=======================================================================
Upon each execution, *Snakemake* automatically creates a log file where all the standard output is recorded. These can be found from the *working directory* into::

    .snakemake/log/

RSP4ABM logs
=======================================================================
In addition to the default *Snakemake*'s logs, *RSP4ABM* create a log directory upon each execution in ::

    logs/<year>/<months>/<day>/<time>/

This directory contains:

- a copy of the executed *Snakemake* command (*cmd.txt*)
- the git commit hash which indicates the version of the RST4ABM (*git.txt*)
- the ID of the user who run the pipeline (*user.txt*)
- a copy of the sample sheet (*local_samples.tsv* or *sra_samples.tsv*)
- a copy of the *config file* (*config.yaml*)

In addition, almost all rules of RST4ABM generate a log file upon execution which records the output of the executed tools or script. These log files are organized in subdirectories of the log directory, mirroring the structure of the main pipeline.  


************************************************************************
Sequencing reads QC
************************************************************************

 :github:`QC rules <rules/0_preprocessing/QC_raw_reads.rules>` assess the sequencing quality of all each sample with FastQC [1]_. Then, a MultiQC [2]_ report generates a report for each sequencing run (based on values of the *sample sheet* column indicated by the "run_column" parameter of the *config file*). A global MultiQC report is generated as well, but without interactive features to deal with the high number of samples  


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





.. _`Snakemake`: https://github.com/metagenlab/microbiome16S_pipeline/tree/master/rules
.. _`Conda environments`: https://github.com/metagenlab/microbiome16S_pipeline/tree/master/envs
.. _`Singularity containers`: https://github.com/metagenlab/microbiome16S_pipeline/tree/master/envs/singularity
.. _`VSEARCH`: https://github.com/torognes/vsearch/releases 
