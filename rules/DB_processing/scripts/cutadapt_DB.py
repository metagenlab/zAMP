from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from snakemake import shell


forward_primer_compl = Seq.reverse_complement(Seq(snakemake.params['forward_primer'], IUPAC.ambiguous_dna))
reverse_primer_compl = Seq.reverse_complement(Seq(snakemake.params['reverse_primer'], IUPAC.ambiguous_dna))
f_length = len(snakemake.params['forward_primer'])
print(snakemake.params['forward_primer'])
print(f_length)
r_length = len(snakemake.params['reverse_primer'])
print(r_length)
print(reverse_primer_compl)

shell("""cutadapt \
--cores {snakemake.threads} \
--error-rate 2 \
--times 1 \
-o {snakemake.output[0]} \
-g '{snakemake.params[forward_primer]}...{reverse_primer_compl}' \
--discard-untrimmed \
{snakemake.input[0]}""")
