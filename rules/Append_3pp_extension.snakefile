# The extensions are coverted into fasta seqs and appended to its original fasta sequence

rule append_three_prime_extension:
    input:
        expand("{three_prime_extensions}/Extended_three_prime_extensions.fasta",three_prime_extensions=three_prime_extensions)

rule do_append_three_prime_extensio:
    input:
        expand("{three_prime_extensions}/Cleaned_three_prime_extensions.bed",three_prime_extensions=three_prime_extensions)
    output:
        expand("{three_prime_extensions}/Extended_three_prime_extensions.fasta",three_prime_extensions=three_prime_extensions)
    threads:
        20
    shell:
        """bash bash_scripts/three_prime/3pp_extension_appender_v2.sh <(sort -u -t$'\\t' -k 6,6 {input} | cut -f1,5,6,9) > {output}   """
