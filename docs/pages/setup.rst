
.. _setup:

########################################################################
Setup and system requirements
########################################################################


.. Note:: Provided command-line examples are given as example for a standard unix terminal in bash.

************************************************************************
Operating system and system resource 
************************************************************************

Operating system
=======================================================================
RST4ABM was designed on *Ubuntu 18.04** (Linux) but should compatible with all system capable of installing the dependencies listed on this page.

RAM memory
=======================================================================
Some tools embedded in RST4ABM can be quite demanding on RAM memory. The actual requirement depends on your dataset and is influenced parameters indicated in the config file. The bottleneck usually is the taxonomic assignment. Factors which can raise your RAM requirements are:
- the number of samples
- the bacterial diversity within your samples
- the number of cores (start by reducing the number of core in snakemake if facing issues)

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

    *Git* is available by default in operating systems. 

    *Test if it is installed*::

        #To test if git is installed, make it print its version. will fail if not installed
        git --version

    If it is not installed, `follow the indications on the git installation page. <https://git-scm.com/downloads>`_.



*Conda*
=======================================================================

    What for?
    -----------------------------------------------------------------------

    `Conda is a convenient python-based package and environment manager. <https://docs.conda.io/en/latest>`_
    It  enables the easy installation of *Snakemake*. Furthermore, it can be used (as an alternative to *Singularity containers*) by *Snakemake* and RSP4ABM to install all the packages required to the execution of the pipeline.


    Install
    -----------------------------------------------------------------------
    `Follow the Miniconda3 installation recommendation <https://docs.conda.io/en/latest/miniconda.html>`_



*Snakemake*
=======================================================================

    What for?
    -----------------------------------------------------------------------
    *RSP4ABM* is a *Snakemake* pipeline. Then, is must be installed and available for each execution of the pipeline. 


    Install
    -----------------------------------------------------------------------
    Follow indications on Snakemake `installation page <https://snakemake.readthedocs.io/en/stable/getting_started/installation.html>`_



************************************************************************
Reference database
************************************************************************

As the last step of setup and before the first execution of the pipeline, a dedicated workflow must be executed to prepare and format the reference taxonomy database. For this, refer to :ref:`DB_preprocessing`. 


