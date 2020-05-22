rule mRNA_filter_five_prime_extension:
    input:
        expand("{five_prime_extensions}/mRNA_filtered_five_prime_extensions.bed",five_prime_extensions=five_prime_extensions),
        expand("{five_prime_extensions}/mRNA_filtered_five_prime_extensions_trans_orientation.bed",five_prime_extensions=five_prime_extensions)


rule do_mRNA_filter_five_prime_extension:
    input:
        expand("{five_prime_extensions}/mRNA_five_prime_extensions.bed",
                                                        five_prime_extensions=five_prime_extensions)
    output:
        expand("{five_prime_extensions}/mRNA_filtered_five_prime_extensions.bed",
                                                        five_prime_extensions=five_prime_extensions)
    threads:
        20
    shell:
        """bash bash_scripts/five_prime/mRNA_5pp_extensions_filter.sh -b {input} > {output}"""

rule do_mRNA_filter_five_prime_extension_trans_orientation:
    input:
        expand("{five_prime_extensions}/mRNA_filtered_five_prime_extensions.bed",
                                                        five_prime_extensions=five_prime_extensions)
    output:
        expand("{five_prime_extensions}/mRNA_filtered_five_prime_extensions_trans_orientation.bed",
                                                        five_prime_extensions=five_prime_extensions)
    threads:
        20
    shell:
        """sort -u -t$'\t' -k 6,6 {input}| cut -f1,5,6,9 > {output}"""
