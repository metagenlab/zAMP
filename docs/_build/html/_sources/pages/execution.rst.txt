
########################################################################
Pipeline execution
########################################################################

.. Note:: Provided command-line examples are given as examples and are valid for a standard unix bash terminal.

.. Note:: The instructions hereunder assume that the pipeline setup was completed (see :ref:`setup`), including the preprocessing of the taxonomy reference database (see :ref:`DB_preprocessing`).


************************************************************************
1. Select a method to define sequencing read files as input
************************************************************************

Before executing RSP4ABM, you need to define how you provide the sequencing read files (in *.fastq.gz* format) to the pipeline. The pipeline accepts four methods to define the input read files, with different requirements. Please read hereunder to define the methods best fitting your needs. 

A. Sample column
=======================================================================

    Rational:
    -----------------------------------------------------------------------
    The easy, but not completely foolproof, way to run the pipeline. 

    Working principle:
    -----------------------------------------------------------------------
    A `regex <https://en.wikipedia.org/wiki/Regular_expression>`_-based script automatically matches the sample names with the sequencing read files format located in a directory. 
        
    Requirements:
    -----------------------------------------------------------------------
    - The leftmost column of the *sample sheet* is named "*Sample*" and contains the sample identifiers.  

    - All samples have unique, simple names starting by a letter (and not a number).

        *e.g. "1" will fail as sample name*

    - All samples names are contained without ambiguity within the filename of the *.fastq.gz* read files.

        *e.g. the pipeline will be confused by "Sample2" and "Sample2_bis"*.     

    - All reads are located in one single directory.

        *reads can be actually stored in the directory or be represented by symbolic links (recommended)*

    - The path of the directory containing the reads is designated by the "*link_directory*" parameter in the *config file* ("*links/*" by default). 

    - Reads are considered to be paired-end (the pipeline looks for "*R1*" and "*R2*" *.fastq.gz* files) by default. For single-end or mixed reads, the *sample sheet* must contain a "*LibraryLayout*" column indicating "*single*" (or "*paired*").

    - The path to the *sample sheet* describing your samples is indicated by the "*local_samples*" parameter in the *config file*.


B. OldSampleName column
=======================================================================

    Rational:
    -----------------------------------------------------------------------
    The sample names used for sequencing and contained in the *.fastq.gz* read filenames can be impractical or incompatible to be used as sample identifier (e.g. stars with a number). In this case, match can be done with an "*OldSampleName*" column instead of the default "*Sample*" column. 

    Working principle:
    -----------------------------------------------------------------------
    Same than for `A. Sample column`_ matching except that the presence of a column named "*OldSampleName*" in the *sample sheet* forces the pipeline to run the `regex <https://en.wikipedia.org/wiki/Regular_expression>`_-based matching on this column instead of the leftmost "*Sample*" column. 

    This method implies that sample are given a "new" identifier from the *Sample* leftmost column of the *sample sheet*, instead of their "old" identifier (the one in the *fastq.gz* filename and in the *OldSampleName* column).


    Requirements:
    -----------------------------------------------------------------------
    - A column named "*OldSampleName*" where all samples have a unique name contained without ambiguity within their *.fastq.gz* filename.
    - Similar to the *Sample column matching* method for the rest.    


C. Absolute path
=======================================================================

    Rational:
    -----------------------------------------------------------------------
    The foolproof approach, pointing by their absolute path the sequencing read files of each sample. 


    Working principle:
    -----------------------------------------------------------------------
    In presence of a "*R1*" column in the *sample sheet* (+/- a "*R2*" column for paired-end reads), the pipeline considers the path indicated in this/these column(s) as the input file(s) for each sample. 

    Requirements:
    -----------------------------------------------------------------------
    - A column named "*R1*" (and a "*R2*" for paired-end reads) with the path to the *.fastq.gz* read files for each sample. 

    - The presence of a value in a "*R2*" column indicates to the pipeline that the sample was sequenced in paired-end.

    - The path to the *sample sheet* describing your samples is indicated by the "*local_samples*" parameter in the *config file*.



D. Sample Read Archive (SRA)
=======================================================================

    Rational:
    -----------------------------------------------------------------------
    RSP4ABM is capable to automatically download and process reads stored on the Sequence Read Archive (`SRA <https://en.wikipedia.org/wiki/Sequence_Read_Archive>`_). 

    Working principle:
    -----------------------------------------------------------------------
    Reads stored on the SRA are designated by their *Run* identifier in the *sample sheet*. The sequencing data are automatically downloaded and then processed by the pipeline.  

    Requirements:
    -----------------------------------------------------------------------
    - The leftmost column of the *sample sheet*  is named "*Run*" and matches "*Run*" identifiers from SRA. 

    - A "*LibraryLayout*" column in the *sample sheet* indicates the sequencing layout ("*single*" or "*paired*").

    - The path to the *sample sheet* describing your samples is indicated by the "*sra_samples*" parameter in the *config file*.



************************************************************************
2. Create a working directory
************************************************************************

Before the execution of the pipeline, prepare a new dedicated directory somewhere on the system where you have sufficient space (NOT within the pipeline directory).

*for instance*::

    #make a new directory named "new_analysis"
    $ mkdir new_analysis

    # Move to the new directory
    $ cd new_analysis


The *config file*, *the sample sheet* and eventually the *links/* directory (see below) must be created within this working directory


******************************************************************************
3. Create a *links* directory (*for Sample and OldSampleName column as input*)
******************************************************************************
When using the method of the `A. Sample column`_ or the `B. OldSampleName column`_ for the definition of the sequencing read files used as input, all the *.fastq.gz* files must be located within one single directory. By default, this directory is named "*links/*"and is found in the working directory, but another path can be defined by the "*link_directory*" parameter of the *config file*. It is recommended to create symbolic links of the original reads into this *links/* directory.

*for instance*::

    #make a new directory named "links"
    $ mkdir links

    # Create symbolic links for .fastq.gz files from distant directory into "links/"
    $ ln -s path/to/a/distant/directory/containing/some/reads/*.fastq.gz links/


************************************************************************
4. Write a *config file*
************************************************************************
Parameters must be provided to adapt the pipeline to the specificities of the analysis. This is done through a *config file* in the `yaml <https://en.wikipedia.org/wiki/YAML>` format. Create this file in the working directory, copy the content of the example below and adapt it to the analysis. The three first parameters (*"link_directory", *"sra_samples"* and *"local_samples"*) are associated to the definition of input files. For the meaning of the other parameters, please refer to the short comment next by each parameter and the detailed description of the pipeline on the :ref:`under_the_hood` page.

*for instance*::

    # Open a graphic text editor and create a config file. Once opened, copy the example below and adapt it. 
    $ gedit config.yaml

    # Or use a command-line text editor, e.g. 
    # $ nano config.yaml 


**Config file example:**

.. literalinclude:: ../../ressources/template_files/config.yaml
    :language: yaml

************************************************************************
5. Write a *sample sheet* (aka *metadata table*):
************************************************************************
The pipeline requires a spreadsheet in `tabulation-separated values (tsv) format <https://en.wikipedia.org/wiki/Tab-separated_values>`_ listing all the file in the analysis and where:

- The leftmost column ("*Sample*" or "*Run*") is the sample identifier. 

    *This identifier must be unique, start with a letter and not a number and cannot contain spaces or "-". "_" are OK.*  
    
- A "*run_column*" describes the sequencing run of each sample.

    *Each sample must have a different value under this column for each sequencing run included in the analysis. If all samples were sequenced together, then the same value must be repeated for all samples. Prefer an alphanumeric factor, e.g. "run_20200101"*

- A "*grouping_column*" regroups the samples for visualization purpose. 

    *Some of the visualization generated by the pipeline will be generated individually for each value contained in the "grouping column"*

- A "*sample_label*" describes each sample.

    *This column must provide a unique, explicit description of each sample. It can be a replication of the *Sample* column but also provides the opportunity to have more concise or explicit description of each sample*.

- Optional (but recommended) columns describes technical metadata. 

    *In this recommended to provide technical metadata (e.g. library preparation DNA yield) to support technical QC of the data*

- Optional (but recommended) columns describes experimental or clinical metadata.

    *In this recommended to provide clinical or experimental description of each sample, which will support later interpretation of the data.*


*For instance:*

.. literalinclude:: ../../ressources/template_files/example_local_samples.tsv 
    :language: tsv
    :start-after: #Example:


************************************************************************
6. Make sure Snakemake is available
************************************************************************
The recommended way to install *Snakemake* is to create a dedicated `Conda environment <https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>`_. (see :ref:`Setup`). Thus, make sure to activate this environment.

*for instance*::

    # Activate your conda environment
    $ conda activate snakemake

    # Test that Snakemake is available by printing its version 
    $ snakemake --version 



************************************************************************
7. Run the pipeline
************************************************************************

.. Hint:: At this stage, the content of your working directory depends upon what was selected in `1. Select a method to define sequencing read files as input`_. 


*Example of directory content for the basic Sample column strategy*::

    # List the content of your working directory
    $ ls
    config.yaml  links  local_sample.tsv


Then, the execution of the pipeline follows the principles of any `*Snakemake pipeline execution* <https://snakemake.readthedocs.io/en/v5.14.0/executing/cli.html>`_. 

But to make it short, here are the requirement arguments.  

.. literalinclude:: ../../ressources/template_files/basic_snakemake_bash_command.sh 
    :language: bash

