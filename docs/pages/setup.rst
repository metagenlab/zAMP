
.. _setup:

########################################################################
System requirements and setup
########################################################################

************************************************************************
Short description 
************************************************************************

**TO BE WRITTEN**



************************************************************************
Detailed description 
************************************************************************


.. Note:: Provided command-line examples are given as examples and are valid for a standard unix bash terminal.


Operating system and system resource 
=======================================================================

Operating system
-----------------------------------------------------------------------
zAMP was designed on *Ubuntu 18.04* (Linux) but should be compatible with all systems capable of installing the dependencies listed on this page.

RAM memory
-----------------------------------------------------------------------
Some tools embedded in zAMP can be quite demanding on RAM memory. The actual requirement depends on your dataset and is influenced by parameters set in the *config file*. The bottleneck usually is the taxonomic assignment. Factors which can increase the RAM requirement are:

- the number of samples
- the bacterial diversity within your samples
- the number of cores

.. Hint:: In practice, 16 to 32 GB are usually required. 



Software dependencies
=======================================================================

.. _git:

*Git*
-----------------------------------------------------------------------

What for?
_______________________________________________________________________

    *Git* is required to download (clone) zAMP. 


Install
_______________________________________________________________________

    *Git* is available by default in operating systems. If not, `follow the indications on the git installation page. <https://git-scm.com/downloads>`_.

    
Test
_______________________________________________________________________

*To test if git is installed*::

    # To test if git is installed, make it print its version. It will fail if it is not installed
    $ git --version



*Conda*
-----------------------------------------------------------------------

What for?
_______________________________________________________________________

    `Conda is a convenient python-based package and environment manager. <https://docs.conda.io/en/latest>`_
    It enables the easy installation of *Snakemake*. Furthermore, it can be used (as an alternative to *Singularity containers*) by *Snakemake* to retrieve all the packages required for the execution of the zAMP.


Install
_______________________________________________________________________
    `Follow the Miniconda3 installation recommendation <https://docs.conda.io/en/latest/miniconda.html>`_.


Test
_______________________________________________________________________

*To test if Conda is installed*::

    # To test if Conda is installed, make it print its version. It will fail if it is not installed
    conda --version



*Mamba*
-----------------------------------------------------------------------

What for?
_______________________________________________________________________

    `Mamba is an alternative to standard conda managers which  <https://docs.conda.io/en/latest>`_
    It enables the easy installation of *Snakemake*. Furthermore, it can be used (as an alternative to *Singularity containers*) by *Snakemake* to retrieve all the packages required for the execution of the zAMP.


Install
_______________________________________________________________________
    `Follow the Mamba installation recommendation <https://github.com/mamba-org/mamba>`_.


*Installed Mamba with Conda*::
    
    mamba install xtensor-r -c conda-forge



Test
_______________________________________________________________________

*To test if Mamba is installed*::

    # To test if Mamba is installed, make it print its version. It will fail if it is not installed
    mamba --version





.. _snakemake:    
 
*Snakemake*
-----------------------------------------------------------------------

What for?
-----------------------------------------------------------------------
    *zAMP* is a *Snakemake* [1]_ pipeline. Therefore, it must be installed and available for execution of the pipeline. 


Install
_______________________________________________________________________
    Follow indications on the *Snakemake* `installation page <https://snakemake.readthedocs.io/en/stable/getting_started/installation.html>`_. It is good practice to create a dedicated *Conda* environment for *Snakemake*. Even if the pipeline should work with newer versions, it was fully tested with Snakemake version 5.26.1. 
    

*To install Snakemake in a dedicated "Snakemake" environment*::

    # Install Snakemake version 5.26.1 in a environment named "snakemake5261"
    mamba create -c bioconda -n snakemake5261  snakemake=5.26.1


Test
_______________________________________________________________________

*To test if Snakemake is installed*::

    # To test if Snakmeake is installed, make it print its version. It will fail if it is not installed
    snakemake --version


.. _singularity:   

*Singularity* 
-----------------------------------------------------------------------

What for?
-----------------------------------------------------------------------

    *Singularity* is a container plateform. It enables to create, retrieve and install containers, which are predefined transposable sets of software. The installation of *Singularity* is optional for most of the functions in zAMP except for the :ref:`*in silicod* prediction pipeline for which it is a requirement <insilico>`. Indeed, the user can choose either Conda_ or Singularity_ to retrieve all the required tools. Yet, it is recommended running zAMP with *Singularity* containers since it enables the best level of reproducibility [2]_. 

    
Install
_______________________________________________________________________
    Follow indications on *Singularity* `installation page <https://sylabs.io/guides/3.6/user-guide/quick_start.html#quick-installation-steps>`_


Test
_______________________________________________________________________

*To test if Singularity is installed*::

    # To test if Singularity is installed, make it print its version. It will fail if it is not installed
    singularity --version




Clone zAMP
=======================================================================

Once all dependencies are installed and working, zAMP can be cloned with git::

    git clone https://github.com/metagenlab/microbiome16S_pipeline.git --recursive


.. Hint:: Please note the path of the directory in which you cloned zAMP since you will need it to execute the pipeline. 



Reference database
=======================================================================

In the very last step of setup and before the first execution of the pipeline, a dedicated workflow must be executed to prepare and format the reference taxonomy database. For this, refer to :ref:`DB_preprocessing` page. 



References
=======================================================================

.. [1] Köster J, Rahmann S. Snakemake-a scalable bioinformatics workflow engine. Bioinformatics. 2012. 
.. [2] Grüning B, Chilton J, Köster J, Dale R, Soranzo N, van den Beek M, et al. Practical Computational Reproducibility in the Life Sciences. Cell Systems. 2018. 
