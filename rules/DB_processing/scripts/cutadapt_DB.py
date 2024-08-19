## Import modules
from Bio.Seq import Seq
from snakemake import shell

## Compute some values to be fed to cutadapt
if snakemake.params['amplicon_type'] == "ITS":
    forward_primer_compl = Seq.reverse_complement(Seq(snakemake.params['forward_primer']))
    reverse_primer_compl = Seq.reverse_complement(Seq(snakemake.params['reverse_primer']))
    f_length = len(snakemake.params['forward_primer'])
    r_length = len(snakemake.params['reverse_primer'])
    f_min_coverage = round(f_length*snakemake.params['coverage'])
    r_min_coverage = round(r_length*snakemake.params['coverage'])

    ## Run cutadapt in shell
    shell("""cutadapt \
    --cores {snakemake.threads} \
    --error-rate {snakemake.params[excepted_errors]} \
    --times 1 \
    -o {snakemake.output[0]} \
    -g '{snakemake.params[forward_primer]};min_overlap={f_min_coverage};optional...{reverse_primer_compl};min_overlap={r_min_coverage}' \
    --discard-untrimmed \
    --minimum-length {snakemake.params[min_length]} \
    --maximum-length {snakemake.params[max_length]} \
    {snakemake.input[0]} >> {snakemake.log[0]}""")

elif snakemake.params['amplicon_type'] == "16S":
    forward_primer_compl = Seq.reverse_complement(Seq(snakemake.params['forward_primer']))
    reverse_primer_compl = Seq.reverse_complement(Seq(snakemake.params['reverse_primer']))
    f_length = len(snakemake.params['forward_primer'])
    r_length = len(snakemake.params['reverse_primer'])
    f_min_coverage = round(f_length*snakemake.params['coverage'])
    r_min_coverage = round(r_length*snakemake.params['coverage'])

    ## Run cutadapt in shell
    shell("""cutadapt \
    --cores {snakemake.threads} \
    --error-rate {snakemake.params[excepted_errors]} \
    --times 1 \
    -o {snakemake.output[0]} \
    -g '{snakemake.params[forward_primer]};min_overlap={f_min_coverage}...{reverse_primer_compl};min_overlap={r_min_coverage}' \
    --discard-untrimmed \
    --minimum-length {snakemake.params[min_length]} \
    --maximum-length {snakemake.params[max_length]} \
    {snakemake.input[0]} >> {snakemake.log[0]}""")
