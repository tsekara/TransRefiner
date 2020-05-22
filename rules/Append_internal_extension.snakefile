# The common 3pp and 5pp extensions are converted into fasta seqs and appended to its original fasta sequence

rule append_internal_extension:
    input:
        expand("{internal_extensions}/Extended_internal_extensions.fasta",internal_extensions=internal_extensions)

rule do_append_internal_extension:
    input:
        five_pp_extn=expand("{internal_extensions}/Cleaned_five_prime_extensions.bed",internal_extensions=internal_extensions),
        three_pp_extn=expand("{internal_extensions}/Cleaned_three_prime_extensions.bed",internal_extensions=internal_extensions)
    output:
        expand("{internal_extensions}/Extended_internal_extensions.fasta",internal_extensions=internal_extensions)
    threads:
        20
    shell:
        """bash bash_scripts/internal/internal_extensions_appender.sh <(awk 'NR==FNR {{ seen[$6]++; next }} seen[$6]' {input.five_pp_extn} {input.three_pp_extn} | cut -f1,5,6,9 | sort -u -t$'\\t' -k3,3) > {output}   """
