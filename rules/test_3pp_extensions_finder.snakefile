rule test_three_prime_extension:
    input:
        expand("{three_prime_extensions}/test_three_prime_extensions.bed",
                                        three_prime_extensions=three_prime_extensions,
                                        RPF_bam=Tophat2_RPF_bam,
                                        search_length=search_length)

rule do_test_three_prime_extension:
    input:
        trans_loc=expand("{transcript_location_in_genome}/3prime_partial_transcripts_in_genome.txt",
                                        transcript_location_in_genome=transcript_location_in_genome),

        annotation=expand("{Gmap_dir}/GMAP_3prime_partial.{extn}",
                                                        extn="gff3",
                                                        Gmap_dir=Gmap_dir)
    output:
        expand("{three_prime_extensions}/test_three_prime_extensions.bed",
                                three_prime_extensions=three_prime_extensions)
    params:
        RPF_bam=Tophat2_RPF_bam,
        search_length=search_length
    threads:
        20
    shell:
        """bash bash_scripts/three_prime/test_3pp_extensions_finder.sh \
        --bam {params.RPF_bam} \
        --search_length {params.search_length} \
        --annotation {input.annotation} \
        --orientation {input.trans_loc} > {output}   """
