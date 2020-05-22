rule test_append_three_prime_extension:
    input:
        expand("{three_prime_extensions}/test_Extended_three_prime_extensions.fasta",
        three_prime_extensions=three_prime_extensions,
        genome_fasta=genome,
        genome_dir=genome_dir,
        transcriptome_fasta=transcriptome)

rule do_test_append_three_prime_extension:
    input:
        filtered_extensions=expand("{three_prime_extensions}/test_filtered_three_prime_extensions.bed",
                                    three_prime_extensions=three_prime_extensions),
        transcript_orientation=expand("{three_prime_extensions}/test_filtered_three_prime_extensions_trans_orientation.bed",
                                    three_prime_extensions=three_prime_extensions)
    output:
        expand("{three_prime_extensions}/test_Extended_three_prime_extensions.fasta",
                                    three_prime_extensions=three_prime_extensions)
    threads:
        20
    params:
        genome_fasta=genome,
        transcriptome_fasta=transcriptome
    shell:
        """bash bash_scripts/three_prime/test_3pp_extensions_appender.sh \
         -b {input.filtered_extensions} \
         -gfa {params.genome_fasta} \
         -tfa {params.transcriptome_fasta} \
         -o {input.transcript_orientation} > {output}"""
