rule Clean_five_prime_extension:
    input:
        expand("{five_prime_extensions}/Cleaned_five_prime_extensions.bed",five_prime_extensions=five_prime_extensions)

rule do_Clean_five_prime_extension:
    input:
        expand("{five_prime_extensions}/five_prime_extensions.bed",five_prime_extensions=five_prime_extensions)
    output:
        expand("{five_prime_extensions}/Cleaned_five_prime_extensions.bed",five_prime_extensions=five_prime_extensions)
    threads:
        20
    shell:
        """bash bash_scripts/five_prime/checking_counts_and_strand_orientation.sh {input} > {output}   """
