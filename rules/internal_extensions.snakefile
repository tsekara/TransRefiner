rule internal_extensions:
    input:
        expand("{internal_extensions}/three_prime_extensions.bed",internal_extensions=internal_extensions),
        expand("{internal_extensions}/five_prime_extensions.bed",internal_extensions=internal_extensions)

rule do_internal_three_prime_extensions:
    input:
        expand("{transcript_location_in_genome}/internal_transcripts_in_genome.txt",transcript_location_in_genome=transcript_location_in_genome)
    output:
        expand("{internal_extensions}/three_prime_extensions.bed",internal_extensions=internal_extensions)
    threads:
        20
    shell:
        """bash bash_scripts/internal/3pp_RPF_extensions_finder_intron_length_200.sh {input} > {output} """

rule do_internal_five_prime_extensions:
    input:
        expand("{transcript_location_in_genome}/internal_transcripts_in_genome.txt",transcript_location_in_genome=transcript_location_in_genome)
    output:
        expand("{internal_extensions}/five_prime_extensions.bed",internal_extensions=internal_extensions)
    threads:
        20
    shell:
        """bash bash_scripts/internal/5pp_RPF_extensions_finder_intron_length_200.sh {input} > {output}"""
