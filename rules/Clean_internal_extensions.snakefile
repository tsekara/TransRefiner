rule internal_Cleaned_extensions:
    input:
        expand("{internal_extensions}/Cleaned_three_prime_extensions.bed",internal_extensions=internal_extensions),
        expand("{internal_extensions}/Cleaned_five_prime_extensions.bed",internal_extensions=internal_extensions)

rule do_Clean_internal_three_prime_extension:
    input:
        expand("{internal_extensions}/three_prime_extensions.bed",internal_extensions=internal_extensions)
    output:
        expand("{internal_extensions}/Cleaned_three_prime_extensions.bed",internal_extensions=internal_extensions)
    threads:
        20
    shell:
        """bash bash_scripts/three_prime/3pp_checking_counts_and_strand_orientation.sh {input} > {output}   """

rule do_Clean_internal_five_prime_extension:
    input:
        expand("{internal_extensions}/five_prime_extensions.bed",internal_extensions=internal_extensions)
    output:
        expand("{internal_extensions}/Cleaned_five_prime_extensions.bed",internal_extensions=internal_extensions)
    threads:
        20
    shell:
        """bash bash_scripts/five_prime/checking_counts_and_strand_orientation.sh {input} > {output}   """
