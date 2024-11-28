Frequently asked questions (FAQ)
=======================================

##1.General Questions

###Who should use zAMP?

Researchers working on microbiota profiling using 16S or ITS amplicon sequencing can use zAMP to process, analyze, and visualize their data efficiently.

##2.Metadata Issues

###Why does zAMP fail to read my metadata file?

Ensure the metadata file adheres to the required format:

- **File format**: TSV (tab-separated values).

- **Mandatory columns**: Include `seq_run` (identifying sequencing runs) and `Sample` (unique identifiers for each sample).

- **Avoid space or special characters**: Space or any special character in sample names cause errors.
