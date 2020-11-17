from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from snakemake import shell


forward_primer_compl = Seq.reverse_complement(Seq(snakemake.params['forward_primer'], IUPAC.ambiguous_dna))
reverse_primer_compl = Seq.reverse_complement(Seq(snakemake.params['reverse_primer'], IUPAC.ambiguous_dna))


shell("""cutadapt \
--cores {snakemake.threads} \
--error-rate 0.3 \
--overlap 5 \
-l {snakemake.params[length]}} \
-o {snakemake.output[0]} \
-g {snakemake.params[reverse_primer]}...{reverse_primer_compl} \
--discard-untrimmed \
{snakemake.input[0]} >> {snakemake.log[0]}""")
