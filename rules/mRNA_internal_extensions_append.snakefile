rule mRNA_append_internal_extension:
    input:
        expand("{internal_extensions}/mRNA_Extended_internal.fasta",internal_extensions=internal_extensions)

rule do_mRNA_append_internal_extension:
    input:
        filtered_5pp_extns=expand("{internal_extensions}/mRNA_filtered_internal_5pp_extensions.bed",
                                                                        internal_extensions=internal_extensions),
        filtered_3pp_extns=expand("{internal_extensions}/mRNA_filtered_internal_3pp_extensions.bed",
                                                                        internal_extensions=internal_extensions),
        internal_trans_loc=expand("{internal_extensions}/mRNA_filtered_internal_5pp_3pp_extensions_trans_orientation.bed",
                                                                        internal_extensions=internal_extensions)
    output:
        expand("{internal_extensions}/mRNA_Extended_internal.fasta",internal_extensions=internal_extensions)
    threads:
        20
    params:
        genome_fasta=genome,
        transcriptome_fasta=transcriptome
    shell:
        """bash bash_scripts/internal/mRNA_internal_extensions_appender.sh \
         -b_5pp {input.filtered_5pp_extns} \
         -b_3pp {input.filtered_3pp_extns} \
         -gfa {params.genome_fasta} \
         -tfa {params.transcriptome_fasta} \
         -o {input.internal_trans_loc} > {output}"""
