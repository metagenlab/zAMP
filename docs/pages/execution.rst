
.. _execution:

########################################################################
Running zAMP
########################################################################

.. Note:: Before running zAMP see :ref:`setup` and :ref:`tax_DB` sections.


zAMP can use reads locally or from NCBI's `Sequence Read Archive (SRA) <https://www.ncbi.nlm.nih.gov/sra>`_.

Local reads
-----------

Sample sheet only
~~~~~~~~~~~~~

*Command*::

    zamp run -i samples.tsv


*Example samples.tsv* :

.. csv-table:: 
    :header-rows: 1

    sample,R1,R2,sample_group,run
    SRR9067116,reads/SRR9067116_1.fastq.gz,reads/SRR9067116_2.fastq.gz,Genital_tract,run1
    SRR9067115,reads/SRR9067115_1.fastq.gz,reads/SRR9067115_2.fastq.gz,Genital_tract,run1
    SRR9067114,reads/SRR9067114_1.fastq.gz,reads/SRR9067114_2.fastq.gz,Genital_tract,run1
    SRR7225909,reads/SRR7225909_1.fastq.gz,reads/SRR7225909_2.fastq.gz,human_biliary_tract,run2
    SRR7225908,reads/SRR7225908_1.fastq.gz,reads/SRR7225908_2.fastq.gz,human_biliary_tract,run2
    SRR7225907,reads/SRR7225907_1.fastq.gz,reads/SRR7225907_2.fastq.gz,human_biliary_tract,run2

* `sample`: the name of the sample
* `R1`: path to forward reads
* `R2`: path to reverse reads
* `sample_group`: sample grouping for visualizations
* `run`: column for applying DADA2 error learning and denoising for each sequencing run 


.. Note:: You can add any other columns in the table provided the above mentionned columns are present

Reads directory and metadata as input
~~~~~~~~~~~~~
*Command*::

    zamp run -i reads -m metadata.tsv


*Example metadata.tsv* :

.. csv-table::
    :header-rows: 1

    fastq,sample,sample_group,run
    SRR9067116,Vaginal-Library42,Genital_tract,run1
    SRR9067115,Vaginal-Library41,Genital_tract,run1
    SRR9067114,Vaginal-Library48,Genital_tract,run1
    SRR7225909,NE14,human_biliary_tract,run2
    SRR7225908,A3D12,human_biliary_tract,run2
    SRR7225907,NN15,human_biliary_tract,run2



SRA reads
-----------

zAMP can fetch reads from NCBI's `SRA <https://www.ncbi.nlm.nih.gov/sra>`_ using `fasterq-dump <https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump>`_.

*Command*::

    zamp run -i sra-samples.tsv

*Example sra-samples.tsv* :

.. csv-table::
    :header-rows: 1

    sample,sample_label,sample_group,run,paired
    SRR9067116,Vaginal-16s-V3V4-Library42,Genital_tract,run1,True
    SRR9067115,Vaginal-16s-V3V4-Library41,Genital_tract,run1,True

* `sample_label`: label to rename the SRA fastq files
* `paired`: whether reads are paired or not (required because snakemake can't guess the pairing from the fastq outputs easily )