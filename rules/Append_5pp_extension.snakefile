# The 5pp extensions are coverted into fasta seqs and appended to its original fasta sequence

rule append_five_prime_extension:
    input:
        expand("{five_prime_extensions}/Extended_five_prime_extensions.fasta",five_prime_extensions=five_prime_extensions)

rule do_append_five_prime_extension:
    input:
        expand("{five_prime_extensions}/Cleaned_five_prime_extensions.bed",five_prime_extensions=five_prime_extensions)
    output:
        expand("{five_prime_extensions}/Extended_five_prime_extensions.fasta",five_prime_extensions=five_prime_extensions)
    threads:
        20
    shell:
        """bash bash_scripts/five_prime/5pp_extensions_appender_v3.sh <(sort -u -t$'\\t' -k 6,6 {input} | cut -f1,5,6,9) > {output}   """
