from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from snakemake import shell


forward_primer_compl = Seq.reverse_complement(Seq(snakemake.params['forward_primer'], IUPAC.ambiguous_dna))
reverse_primer_compl = Seq.reverse_complement(Seq(snakemake.params['reverse_primer'], IUPAC.ambiguous_dna))
f_length = len(snakemake.params['forward_primer'])
print(f_length)
r_length = len(snakemake.params['reverse_primer'])
print(r_length)

shell("""cutadapt \
--cores {snakemake.threads} \
--error-rate 0.2 \
--times 1 \
-l {snakemake.params[length]} \
-o {snakemake.output[0]} \
-g '{snakemake.params[forward_primer]};min_overlap=16...{reverse_primer_compl};min_overlap=20' \
--discard-untrimmed \
{snakemake.input[0]} >> {snakemake.log[0]}""")
