from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

    if params.amplicon_type == 'ITS' 
        forward_primer_compl = Seq.reverse_complement(Seq(params.forward_primer, IUPAC.ambiguous_dna)),
        reverse_primer_compl = Seq.reverse_complement(Seq(params.reverse_primer, IUPAC.ambiguous_dna))
        shell("
            cutadapt \
            --cores {threads} \
            --error-rate 0.1 \
            --times 1 \
            --overlap 3 \
            -o {output[R1_trimmed_reads]} \
            -p {output[R2_trimmed_reads]} \
            -g '{params[forward_primer]}' -a '{reverse_primer_compl}' \
            -G '{params[reverse_primer]}' -A '{forward_primer_compl}' \
            --match-read-wildcards \
            {input[R1_raw_reads]} \
            {input[R2_raw_reads]} >> {log}
        ")

    else
        shell("
            cutadapt \
            --cores {threads} \
            --error-rate 0.1 \
            --times 1 \
            --overlap 3 \
            -o {output[R1_trimmed_reads]} \
            -p {output[R2_trimmed_reads]} \
            -g '{params[forward_primer]}' \
            -G '{params[reverse_primer]}' \
            --match-read-wildcards \
            {input[R1_raw_reads]} \
            {input[R2_raw_reads]} >> {log}
        ")