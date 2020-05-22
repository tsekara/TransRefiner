rule mRNA_three_prime_extension_finder:
    input:
        expand("{three_prime_extensions}/mRNA_three_prime_extensions.bed",
                                        three_prime_extensions=three_prime_extensions,
                                        mRNA_bam=Tophat2_mRNA_bam,
                                        search_length=search_length)

rule do_mRNA_three_prime_extension_finder:
    input:
        trans_loc=expand("{transcript_location_in_genome}/3prime_partial_transcripts_in_genome.txt",
                                        transcript_location_in_genome=transcript_location_in_genome),

        annotation=expand("{Gmap_dir}/GMAP_3prime_partial.{extn}",
                                                        extn="gff3",
                                                        Gmap_dir=Gmap_dir)
    output:
        expand("{three_prime_extensions}/mRNA_three_prime_extensions.bed",
                                three_prime_extensions=three_prime_extensions)
    params:
        mRNA_bam=Tophat2_mRNA_bam,
        search_length=search_length
    threads:
        20
    shell:
        """bash bash_scripts/three_prime/mRNA_3pp_extensions_finder.sh \
        --bam {params.mRNA_bam} \
        --search_length {params.search_length} \
        --annotation {input.annotation} \
        --orientation {input.trans_loc} > {output}   """
