rule mRNA_filter_internal_extension:
    input:
        expand("{internal_extensions}/mRNA_filtered_internal_5pp_3pp_extensions_trans_orientation.bed",
                                                    internal_extensions=internal_extensions)
rule do_mRNA_filter_internal_5pp_extension:
    input:
        expand("{internal_extensions}/mRNA_internal_5pp_extensions.bed",
                                internal_extensions=internal_extensions)
    output:
        expand("{internal_extensions}/mRNA_filtered_internal_5pp_extensions.bed",
                                        internal_extensions=internal_extensions)
    threads:
        20
    shell:
        """bash bash_scripts/five_prime/mRNA_5pp_extensions_filter.sh -b {input} > {output}"""


rule do_mRNA_filter_internal_3pp_extension:
    input:
        expand("{internal_extensions}/mRNA_internal_3pp_extensions.bed",
                                internal_extensions=internal_extensions)
    output:
        expand("{internal_extensions}/mRNA_filtered_internal_3pp_extensions.bed",
                                internal_extensions=internal_extensions)
    threads:
        20
    shell:
        """bash bash_scripts/three_prime/mRNA_3pp_extensions_filter.sh -b {input} > {output}"""

rule do_filter_internal_trans_location:
    input:
        five_prime=expand("{internal_extensions}/mRNA_filtered_internal_5pp_extensions.bed",
                                        internal_extensions=internal_extensions),
        three_prime=expand("{internal_extensions}/mRNA_filtered_internal_3pp_extensions.bed",
                                internal_extensions=internal_extensions)
    output:
        expand("{internal_extensions}/mRNA_filtered_internal_5pp_3pp_extensions_trans_orientation.bed",
                                                internal_extensions=internal_extensions)
    shell:
        """awk 'NR==FNR {{ seen[$6]++; next }} seen[$6]' {input.five_prime} {input.three_prime} | cut -f1,5,6,9 | sort -u -t$'\t' -k3,3 > {output}"""

#rule do_mRNA_filter_internal_5pp_extension_trans_orientation:
#    input:
#        expand("{internal_extensions}/mRNA_filtered_internal_5pp_extensions.bed",
#                                        internal_extensions=internal_extensions)
#    output:
#        expand("{internal_extensions}/mRNA_filtered_internal_5pp_extensions_trans_orientation.bed",
#                                        internal_extensions=internal_extensions)
#    threads:
#        20
#    shell:
#        """sort -u -t$'\t' -k 6,6 {input}| cut -f1,5,6,9 > {output}"""



#rule do_mRNA_filter_internal_3pp_extension_trans_orientation:
#    input:
#        expand("{internal_extensions}/mRNA_filtered_internal_3pp_extensions.bed",
#                                        internal_extensions=internal_extensions)
#    output:
#        expand("{internal_extensions}/mRNA_filtered_internal_3pp_extensions_trans_orientation.bed",
#                                        internal_extensions=internal_extensions)
#    threads:
#        20
#    shell:
#        """sort -u -t$'\t' -k 6,6 {input}| cut -f1,5,6,9 > {output}"""
