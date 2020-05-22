rule mRNA_five_prime_extension_finder:
    input:
        expand("{five_prime_extensions}/mRNA_five_prime_extensions.bed",
                                        five_prime_extensions=five_prime_extensions,
                                        mRNA_bam=Tophat2_mRNA_bam,
                                        search_length=search_length)

rule do_mRNA_five_prime_extension_finder:
    input:
        trans_loc=expand("{transcript_location_in_genome}/5prime_partial_transcripts_in_genome.txt",
                                        transcript_location_in_genome=transcript_location_in_genome),

        annotation=expand("{Gmap_dir}/GMAP_5prime_partial.{extn}",
                                                        extn="gff3",
                                                        Gmap_dir=Gmap_dir)
    output:
        expand("{five_prime_extensions}/mRNA_five_prime_extensions.bed",
                                five_prime_extensions=five_prime_extensions)
    params:
        mRNA_bam=Tophat2_mRNA_bam,
        search_length=search_length
    threads:
        20
    shell:
        """bash bash_scripts/five_prime/mRNA_5pp_extensions_finder.sh \
        --bam {params.mRNA_bam} \
        --search_length {params.search_length} \
        --annotation {input.annotation} \
        --orientation {input.trans_loc} > {output}   """
