Frequently asked questions (FAQ)
=======================================


Who should use zAMP?
--------------------
Researchers working on microbiota profiling using 16S or ITS amplicon sequencing can use zAMP to process, analyze, and visualize their data efficiently.



Why does zAMP fail to read my metadata file?
-------------------------------------------

Ensure the metadata file adheres to the required format:

- **File format**: TSV (tab-separated values).

- **Mandatory columns**: Make sure to have mandatory columns like `run`, `sample` and `sample_group` in your sample sheet (see :ref:`execution`)

- **Avoid space or special characters**: Space or any special character in sample names cause errors.



How do I check the available commands and options in zAMP?
-----------------------------------------------------------
You can use the `zamp -h` command to see the general options and a list of available commands:


What parameters must I modify in the command line before running zAMP?
----------------------------------------------------------------------

Required parameters:

- `-i`:  can be either a sample sheet or the path to your reads

- `-m`: Path to your metadata if your input is a reads directory

- `--fw-primer` and `--rv-primer`: for forward and reverse primer sequences.

- `-db`: Path to your database files

What happens if I set incorrect primer sequences?
-------------------------------------------------

Incorrect primers may cause failure during the trimming step or incorrect merging of paired-end reads. Double-check primer sequences for accuracy.



What should I do if zAMP crashes midway?
----------------------------------------

- Identify the error message in the log files.

- Check file paths and parameter settings.




What reference databases are supported by zAMP?
-----------------------------------------------

zAMP supports multiple databases, including:

- SILVA

- Greengenes2

- EzBioCloud




