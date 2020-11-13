
########################################################################
Pipeline execution
########################################################################

************************************************************************
1. Select a method to define sequencing read files as input
************************************************************************

Before executing RSP4ABM, you need to define how you will provide the sequencing read files (in *.fastq.gz* format) to the pipeline. The pipeline accepts four methods to define the input read files, with different requirements. Please read hereunder to define the methods best fitting your needs. 

A. Sample column
=======================================================================

Rational:
-----------------------------------------------------------------------
    The easy, but not completely foolproof, way to run the pipeline. 

Working principle:
-----------------------------------------------------------------------
    A `regex <https://en.wikipedia.org/wiki/Regular_expression>`_-based script will automatically match the names of the samplesfrom the *sample sheet* with the read files in *.fastq.gz* format found in a directory. 
    
Requirements:
-----------------------------------------------------------------------

    - The leftmost column of your *sample sheet* is named "*Sample*" and contains the sample identifiers.  

    - All samples have unique, simple names starting by a letter (and not a number).
        *e.g. "1" will fail as sample name*

    - All samples names are contained without ambiguity within the filename of your *.fastq.gz* read files.
        *e.g. the pipeline will be confused by "Sample2" and "Sample2_bis"*.     

    - All reads are located in one single directory.
        *reads can be actually stored in the folder or be represented by symbolic links (recommended)*

    - The path of the directory containing your reads is designated by the "*link_directory*" parameter in the *config file* ("*links/*" by default). 

    - Reads are considered to be paired-end (the pipeline looks for "*R1*" and "*R2*" *.fastq.gz* files) by default. For single-end or mixed reads, the *sample sheet* must contain a "*LibraryLayout*" column indicating "*single*" (or "*paired*").


B. OldSampleName column
=======================================================================

Rational:
-----------------------------------------------------------------------
    The sample names used for sequencing and contained in the *.fastq.gz* read filenames can impractical or incompatible with our pipeline (e.g. starting with a number). The inclusion of a column named "*OldSampleName*" in the *sample sheet* will force the pipeline to match the sample reads with this column, instead of the *Sample* leftmost column of the spreadsheet.

Working principle:
-----------------------------------------------------------------------
    Same than for *sample column matching* except that the `regex <https://en.wikipedia.org/wiki/Regular_expression>`_-based matching will be done with the "*OldSampleName*" column of the *sample sheet* instead of the leftmost "*Sample*" column.

Requirements:
-----------------------------------------------------------------------
    - A column named "*OldSampleName*" where all samples have a unique name contained without ambiguity within their *.fastq.gz* filename.
    - Similar to the *Sample column matching* method for the rest.    


C. Absolute path
=======================================================================

Rational:
-----------------------------------------------------------------------
    Eventually requires more time to prepare the *sample sheet*. Yet, it is foolproof since the *.fastq.gz* files are designated by the absolute path. 

Working principle:
-----------------------------------------------------------------------
    In presence of a "*R1*" column in the *sample sheet* (+/- a "*R2*" column for paired-end reads), the pipeline will consider the path indicated in this/these column(s) as the input file(s) for each sample. 

Requirements:
-----------------------------------------------------------------------
    - A column named "*R1*" (and a "*R2*" for paired-end reads) with the corresponding path of the *.fastq.gz* read files for each sample. 
    - The presence of a value in a "*R2*" column indicates to the pipeline that the sample was sequenced in paired-end.


D. Sample Read Archive (SRA) hosted reads.
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




************************************************************************
2. Setup your working directory
************************************************************************

.. *Provided command-line examples are for a standard unix terminal in bash* 


0. Working directory
---------------------------------

Before the execution of the pipeline, prepare a new dedicated directory somewhere on your computer where you have sufficient space (and NOT within the pipeline folder).

*for instance*::

    #make a new directory named "new_analysis"
    $ mkdir new_analysis

    # Move to the new directory
    $ cd new_analysis

Then, your working directory must contain the following elements: 
    - config file
    - sample sheet
    - a link directory (optional)

1. Config file:
---------------------------------

Parameters must be provided to adapt the pipeline to the specificities of your analysis. This is done through a config file in the "*.yaml*" format. Create this file in your working directory, copy the content of the example below and adapt it to your analysis. The three first parameters (*"link_directory", *"sra_samples"* and *"local_samples"*) are associated to the definition of input files and are explained hereunder. For the meaning of the other parameters, please refer to the short comment next by each parameter and the detailed description of the pipeline on the :ref:`under_the_hood` page.

*for instance*::

    # Open a graphic text editor and create a config file. Once opened, copy the example below and adapt it. 
    $ gedit config.yaml

    # Or use a command-line text editor, e.g. 
    # $ nano config.yaml 


**Config file example:**

.. literalinclude:: ../../ressources/template_files/config.yaml
    :language: yaml


2. Sample sheet (metadata):
---------------------------------
The pipeline requires a spreadsheet in `tabulation-separated values (tsv) format <https://en.wikipedia.org/wiki/Tab-separated_values>`_ listing all the file in your analysis and where:

- A leftmost "*Sample*" column is an unique sample identifier. 
    *This identifier must be unique, start with a letter and not a number and cannot contain spaces or "-". "_" are OK*  
    
- a "*run_column*" describes the sequencing run of each sample.
    *Each sample must have a different value under this column for each sequencing run included in the analysis. If all samples were sequenced together, then the same value must be repeated for all samples. Prefer an alphanumeric factor, e.g. "run_20200101"*

- a "*grouping_column*" regroups the samples for visualization purpose. 
    *Some of the visualization generated by the pipeline will be generated individually for each value contained in the "grouping column"*

- a "*sample_label*" describes each sample.
    *This column must provide a unique, explicit description of each sample. It can be a replication of the *Sample* column but also provides the opportunity to have more concise or explicit description of each sample*.

- optional (but recommended) columns describes technical metadata. 
    *In this recommended to provide technical metadata (e.g. library preparation DNA yield) to support technical QC of the data*

- optional (but recommended) columns describes experimental or clinical metadata.
    *In this recommended to provide clinical or experimental description of each sample, which will support later interpretation of the data.*


**Sample sheet file example:**

.. literalinclude:: ../../ressources/template_files/example_local_samples.tsv 
    :language: tsv
    :start-after: #Example:



Method 1: local *.fastq* matching 
---------------------------------
With this method, you will deposit all required in a "*links*" folder. The files in this folder can be the actual files or symbolic links to you reads. The pipeline with 
