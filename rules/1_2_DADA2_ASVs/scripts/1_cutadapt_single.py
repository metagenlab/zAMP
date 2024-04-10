from Bio.Seq import Seq
from snakemake import shell

if snakemake.params['amplicon_type'] == "ITS":
    print("ITS Trimming")
    forward_primer_compl = Seq.reverse_complement(Seq(snakemake.params['forward_primer']))
    shell("""cutadapt \
    --cores {snakemake.threads} \
    --error-rate 0.1 \
    --times 2 \
    --overlap 3 \
    -o {snakemake.output[R1_trimmed_reads]} \
    -g '{snakemake.params[forward_primer]}' \
    -a '{reverse_primer_compl}' \
    --match-read-wildcards \
    --discard-untrimmed \
    {snakemake.input[R1_raw_reads]} >> {snakemake.log[0]}""")

elif snakemake.params['amplicon_type'] == "16S":
    print("16S Trimming")
    reverse_primer_compl = Seq.reverse_complement(Seq(snakemake.params['reverse_primer'], IUPAC.ambiguous_dna))
    shell("""cutadapt \
    --cores {snakemake.threads} \
    --error-rate 0.1 \
    --times 1 \
    --overlap 3 \
    -o {snakemake.output[R1_trimmed_reads]} \
    -g '{snakemake.params[forward_primer]}' \
    -a '{reverse_primer_compl}' \
    --match-read-wildcards \
    --discard-untrimmed \
    {snakemake.input[R1_raw_reads]}  >> {snakemake.log[0]}""")




