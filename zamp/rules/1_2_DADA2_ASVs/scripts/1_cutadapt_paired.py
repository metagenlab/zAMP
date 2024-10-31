from Bio.Seq import Seq
from snakemake import shell

print("Primer Trimming")
shell("""cutadapt \
--cores {snakemake.threads} \
--error-rate 0.1 \
--times 1 \
--overlap 3 \
-o {snakemake.output[R1_trimmed_reads]} \
-p {snakemake.output[R2_trimmed_reads]} \
-g '{snakemake.params[forward_primer]}' \
-G '{snakemake.params[reverse_primer]}' \
--match-read-wildcards \
--discard-untrimmed \
{snakemake.input[R1_raw_reads]} \
{snakemake.input[R2_raw_reads]} >> {snakemake.log[0]}""")




