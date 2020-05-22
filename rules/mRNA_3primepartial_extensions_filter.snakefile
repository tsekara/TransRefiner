rule mRNA_filter_three_prime_extension:
    input:
        expand("{three_prime_extensions}/mRNA_filtered_three_prime_extensions.bed",three_prime_extensions=three_prime_extensions),
        expand("{three_prime_extensions}/mRNA_filtered_three_prime_extensions_trans_orientation.bed",three_prime_extensions=three_prime_extensions)


rule do_mRNA_filter_three_prime_extension:
    input:
        expand("{three_prime_extensions}/mRNA_three_prime_extensions.bed",
                                                        three_prime_extensions=three_prime_extensions)
    output:
        expand("{three_prime_extensions}/mRNA_filtered_three_prime_extensions.bed",
                                                        three_prime_extensions=three_prime_extensions)
    threads:
        20
    shell:
        """bash bash_scripts/three_prime/mRNA_3pp_extensions_filter.sh.sh {input} > {output}"""

rule do_mRNA_filter_three_prime_extension_trans_orientation:
    input:
        expand("{three_prime_extensions}/mRNA_filtered_three_prime_extensions.bed",
                                                        three_prime_extensions=three_prime_extensions)
    output:
        expand("{three_prime_extensions}/mRNA_filtered_three_prime_extensions_trans_orientation.bed",
                                                        three_prime_extensions=three_prime_extensions)
    threads:
        20
    shell:
        """sort -u -t$'\t' -k 6,6 {input}| cut -f1,5,6,9 > {output}"""
