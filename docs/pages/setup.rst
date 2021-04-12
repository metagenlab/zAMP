
.. _setup:

########################################################################
Setup and system requirements
########################################################################


.. Note:: Provided command-line examples are given as examples and are valid for a standard unix bash terminal.

************************************************************************
Operating system and system resource 
************************************************************************

Operating system
=======================================================================
RST4ABM was designed on *Ubuntu 18.04* (Linux) but should compatible with all system capable of installing the dependencies listed on this page.

RAM memory
=======================================================================
Some tools embedded in RST4ABM can be quite demanding on RAM memory. The actual requirement depends on your dataset and is influenced by parameters set in the *config file*. The bottleneck usually is the taxonomic assignment. Factors which can increase the RAM requirement are:

- the number of samples
- the bacterial diversity within your samples
- the number of cores

.. Hint:: In practice, 16 to 32 GB are usually required. 


************************************************************************
Software dependencies
************************************************************************

*Git*
=======================================================================

What for?
-----------------------------------------------------------------------

    *Git* is required to download (clone) RSP4ABM. 


Install
-----------------------------------------------------------------------

    *Git* is available by default in operating systems. If not, `follow the indications on the git installation page. <https://git-scm.com/downloads>`_.

    
Test
-----------------------------------------------------------------------

*To test if git is installed*::

    # To test if git is installed, make it print its version. It will fail if it is not installed
    $ git --version



*Conda*
=======================================================================

What for?
-----------------------------------------------------------------------

    `Conda is a convenient python-based package and environment manager. <https://docs.conda.io/en/latest>`_
    It enables the easy installation of *Snakemake*. Furthermore, it can be used (as an alternative to *Singularity containers*) by *Snakemake* to retrieve all the packages required for the execution of the RSP4ABM.


Install
-----------------------------------------------------------------------
    `Follow the Miniconda3 installation recommendation <https://docs.conda.io/en/latest/miniconda.html>`_.


Test
-----------------------------------------------------------------------

*To test if Conda is installed*::

    # To test if Conda is installed, make it print its version. It will fail if it is not installed
    conda --version



*Snakemake*
=======================================================================

What for?
-----------------------------------------------------------------------
    *RSP4ABM* is a *Snakemake* [1]_ pipeline. Therefore, it must be installed and available for execution of the pipeline. 


Install
-----------------------------------------------------------------------
    Follow indications on *Snakemake* `installation page <https://snakemake.readthedocs.io/en/stable/getting_started/installation.html>`_. It is good practice to create a dedicated *Conda* environment for *Snakemake*.
    

*To install Snakemake in a dedicated "Snakemake" environment*::

    # Install snakemake in a environment named "Snakemake"
    conda create -c bioconda -n snakemake  snakemake


Test
-----------------------------------------------------------------------

*To test if Snakemake is installed*::

    # To test if Snakmeake is installed, make it print its version. It will fail if it is not installed
    snakemake --version


*Singularity* 
=======================================================================

What for?
-----------------------------------------------------------------------
    *Singularity* is a container plateform. It enables to create, retrieve and install containers, which are predefined transposable sets of software. The installation of *Singularity* is optional for most of the functions in RSP4ABM. Indeed, the user can choose either Conda_ or Singularity_ to retrieve all the required tools. Yet, it is recommended running RSP4ABM with *Singularity* containers since it enables the best level of reproducibility [2]_. 

    
Install
-----------------------------------------------------------------------
    Follow indications on *Singularity* `installation page <https://sylabs.io/guides/3.6/user-guide/quick_start.html#quick-installation-steps>`_


Test
-----------------------------------------------------------------------

*To test if Singularity is installed*::

    # To test if Singularity is installed, make it print its version. It will fail if it is not installed
    singularity --version



************************************************************************
Clone RSP4ABM
************************************************************************

Once all dependencies installed and working, RSP4ABM can be cloned with git::

    git clone https://github.com/metagenlab/microbiome16S_pipeline.git


Please note the path of the directory in which you cloned RSP4ABM since you will have to indicate it when executing the pipeline. 


************************************************************************
Reference database
************************************************************************

The very last step of setup and before the first execution of the pipeline, a dedicated workflow must be executed to prepare and format the reference taxonomy database. For this, refer to :ref:`DB_preprocessing`. 


************************************************************************
References
************************************************************************

.. [1] Köster J, Rahmann S. Snakemake-a scalable bioinformatics workflow engine. Bioinformatics. 2012. 
.. [2] Grüning B, Chilton J, Köster J, Dale R, Soranzo N, van den Beek M, et al. Practical Computational Reproducibility in the Life Sciences. Cell Systems. 2018. 
