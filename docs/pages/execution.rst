Pipeline execution
==================

*Provided command-line examples a standard unix terminal in bash* 

0. Working directory
---------------------------------

Before the execution of the pipeline, prepare a new dedicated directory

*for instance*::

    #make a new directory named "new_analysis"
    $ mkdir new_analysis

    # change to the new directory
    $ cd new_analysis

Then, your working directory must contain the following elements: 

1. A *config* file:
---------------------------------

Parameters must be provided to adapt the pipeline to the specificities of your analysis. This is done through a config file in the "*.yaml*" format. Create this file in your working directory, copy the content of the example below and adapt it to your analysis. Refer to :ref:`under_the_hood.rst`_ page for detailed meaning of each parameter.


**Config file example:**

.. literalinclude:: ../../ressources/template_files/config.yaml
    :language: yaml




2. A *sample_table*:
---------------------------------
The pipeline requires a spread-sheet in tabulation-separated values (tsv) format listing all the file in your analysis and where:
- the first column in a Unique Sample ID. This ID can NOT start by a number. 
- sufficient column to describe all the parameters demanded in the config file. 
- all metadata


https://github.com/metagenlab/microbiome16S_pipeline/blob/master/ressources/template_files/config.yaml


Method 1: local *.fastq* matching 
---------------------------------
With this method, you will deposit all required in a "*links*" folder. The files in this folder can be the actual files or symbolic links to you reads. The pipeline with 
