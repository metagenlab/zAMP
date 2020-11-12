Setup and Requirements
==========

Operating system and system resource 
----------------------

Operating system
^^^^^^^^^^^^^^^^^^^^^^^^
RST4ABM was designed on Linux 18.04 but should compatible with all system capable of installing the dependencies listed on this page.

RAM memory
^^^^^^^^^^^^^^^^^^^^^^^^
Some tools embedded in RST4ABM can be quite demanding regarding RAM memory. The actual requirements depends on your dataset, are influences by parameters and variates along processing steps. The bottleneck usually is the taxonomic assignment. Factors which can raise your RAM requirements are:
- the number of samples
- the bacterial diversity within your samples
- the number of cores (start by reducing the number of core in snakemake if facing issues)
In practice, 16 to 32 GB are usually required. 

Snakemake
----------
Snakemake is the cornerstone of this pipeline. It can be easily installed through 

For full description of installation options,  look at the dedicated `Snakemake installation page <https://snakemake.readthedocs.io/en/stable/getting_started/installation.html>`_. It can   





Conda
----------
Conda is a python-based package manager enabling to easily install environments containing required tools and software. It 


Singularity 
-------------


Taxonomic reference database
-------------------------------
