
.. _setup:

########################################################################
Installation and resource requirements
########################################################################


Quick start
-----------


* Install with conda via `Miniforge  <https://github.com/conda-forge/miniforge>`_::

    conda install zamp


* Install with docker from `Dockerhub <https://hub.docker.com/r/metagenlab/zamp/>`_::
    
    docker pull metagenlab/zamp:latest




Operating system
-----------------------------------------------------------------------
zAMP is available on `Bioconda <https://bioconda.github.io/>`_ which only support Linux and macOS. 

If you use windows, you can still run zAMP via `WSL <https://learn.microsoft.com/en-us/windows/wsl/install>`_.



Installation methods
-----------------------------------------------------------------------

*From source*
=======================================================================
You can install zAMP from source, by cloning the `repository <https://github.com/metagenlab/zAMP>`_::
    
    git clone https://github.com/metagenlab/zAMP
    pip install -e zAMP/

Dependencies:
    * python >=3.11
    * apptainer
    * conda


*Conda*
=======================================================================

You can install zAMP from `Bioconda <https://bioconda.github.io/>`_ with conda installed from `Miniforge  <https://github.com/conda-forge/miniforge>`_::

    # Install conda from Miniforge
    curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
    bash Miniforge3-$(uname)-$(uname -m).sh

    # Install zamp
    conda install zamp

*Containers*
=======================================================================

You can install zAMP by pulling a container image for maximum reproducibility and no extra dependency installation.

We recommend installing zAMP containers via docker or apptainer from any of the container registries below:

* `Dockerhub <https://hub.docker.com/r/metagenlab/zamp/>`_::

    docker pull metagenlab/zamp:latest

* `ghcr <https://github.com/metagenlab/zAMP/pkgs/container/zamp>`_::

    docker pull ghcr.io/metagenlab/zamp:latest

* `biocontainers <https://quay.io/repository/biocontainers/zamp>`_::

    docker pull quay.io/biocontainers/zamp:1.0.0--pyhdfd78af_1


Resource usage
-----------------------------------------------------------------------
Some steps in zAMP can be quite resource-intensive, requiring more RAM. 

For example, if you want to re-train the RDP classifier on SILVA138.1, you might need around 100GB of RAM. 

In terms of duration, and when using the default zAMP module, the bottleneck is usually the DADA2 denoising step, which can take some time with lots of samples.

Actual resource usage like RAM and CPU-time depends on:

- Number of samples
- Sequencing depth of each sample
- Threads set by the user

.. Note:: Usually, a zAMP run duration is less than 1h and does not require more than 32 GB of RAM.