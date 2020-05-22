rule Clean_three_prime_extension:
    input:
        expand("{three_prime_extensions}/Cleaned_three_prime_extensions.bed",three_prime_extensions=three_prime_extensions)

rule do_Clean_three_prime_extension:
    input:
        expand("{three_prime_extensions}/three_prime_extensions.bed",three_prime_extensions=three_prime_extensions)
    output:
        expand("{three_prime_extensions}/Cleaned_three_prime_extensions.bed",three_prime_extensions=three_prime_extensions)
    threads:
        20
    shell:
        """bash bash_scripts/three_prime/3pp_checking_counts_and_strand_orientation.sh {input} > {output}   """
