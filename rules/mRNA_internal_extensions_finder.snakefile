rule mRNA_internal_extension_finder:
    input:
        expand("{internal_extensions}/mRNA_internal_5pp_extensions.bed",
                                        internal_extensions=internal_extensions,
                                        mRNA_bam=Tophat2_mRNA_bam,
                                        search_length=search_length),
        expand("{internal_extensions}/mRNA_internal_3pp_extensions.bed",
                                        internal_extensions=internal_extensions,
                                        mRNA_bam=Tophat2_mRNA_bam,
                                        search_length=search_length)


rule do_mRNA_internal_5pp_extension_finder:
    input:
        trans_loc=expand("{transcript_location_in_genome}/internal_transcripts_in_genome.txt",
                                        transcript_location_in_genome=transcript_location_in_genome),

        annotation=expand("{Gmap_dir}/GMAP_internal.{extn}",
                                        extn="gff3",Gmap_dir=Gmap_dir)
    output:
        expand("{internal_extensions}/mRNA_internal_5pp_extensions.bed",
                                internal_extensions=internal_extensions)
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

rule do_mRNA_internal_3pp_extension_finder:
    input:
        trans_loc=expand("{transcript_location_in_genome}/internal_transcripts_in_genome.txt",
                                        transcript_location_in_genome=transcript_location_in_genome),

        annotation=expand("{Gmap_dir}/GMAP_internal.{extn}",
                                                        extn="gff3",
                                                        Gmap_dir=Gmap_dir)
    output:
        expand("{internal_extensions}/mRNA_internal_3pp_extensions.bed",
                                internal_extensions=internal_extensions)
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
