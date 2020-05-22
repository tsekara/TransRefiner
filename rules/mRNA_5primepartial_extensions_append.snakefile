rule mRNA_append_five_prime_extension:
    input:
        expand("{five_prime_extensions}/mRNA_Extended_five_prime_extensions.fasta",
        five_prime_extensions=five_prime_extensions,
        genome_fasta=genome,
        transcriptome_fasta=transcriptome)

rule do_mRNA_append_five_prime_extension:
    input:
        filtered_extensions=expand("{five_prime_extensions}/mRNA_filtered_five_prime_extensions.bed",
                                                                        five_prime_extensions=five_prime_extensions),
        transcript_orientation=expand("{five_prime_extensions}/mRNA_filtered_five_prime_extensions_trans_orientation.bed",
                                                                        five_prime_extensions=five_prime_extensions)
    output:
        expand("{five_prime_extensions}/mRNA_Extended_five_prime_extensions.fasta",
                                    five_prime_extensions=five_prime_extensions)
    threads:
        20
    params:
        genome_fasta=genome,
        transcriptome_fasta=transcriptome
    shell:
        """bash bash_scripts/five_prime/mRNA_5pp_extensions_appender.sh \
         -b {input.filtered_extensions} \
         -gfa {params.genome_fasta} \
         -tfa {params.transcriptome_fasta} \
         -o {input.transcript_orientation} > {output}"""
