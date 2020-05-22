rule five_prime_extension:
    input:
        expand("{five_prime_extensions}/five_prime_extensions.bed",five_prime_extensions=five_prime_extensions)

rule do_five_prime_extension:
    input:
        expand("{transcript_location_in_genome}/5prime_partial_transcripts_in_genome.txt",transcript_location_in_genome=transcript_location_in_genome)
    output:
        expand("{five_prime_extensions}/five_prime_extensions.bed",five_prime_extensions=five_prime_extensions)
    threads:
        20
    shell:
        """bash bash_scripts/five_prime/5pp_RPF_extensions_finder_intron_length_200.sh {input} > {output}"""
