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


## 3. Commands and Parameters

### How do I check the available commands and options in zAMP?

You can use the `zamp -h` command to see the general options and a list of available commands:


### What parameters must I modify in the command line before running zAMP?

At minimum, set:

- `local_samples` or `sra_samples`: Define input sample locations.

- `run_column`: Specify the sequencing run column.

- `forward_primer` and `reverse_primer`: Adjust for non-standard primers.

- `tax_DB_path` and `tax_DB_name`: Ensure the correct reference database path and version are specified.

### What happens if I set incorrect primer sequences?

Incorrect primers may cause failure during the trimming step or incorrect merging of paired-end reads. Double-check primer sequences for accuracy.

## 4. Pipeline Execution

### What should I do if zAMP crashes midway?

- Identify the error message in the log files.

- Check file paths and parameter settings.


## 5. Taxonomic Classification

### What reference databases are supported by zAMP?

zAMP supports multiple databases, including:

- SILVA

- Greengenes2

- EzBioCloud




