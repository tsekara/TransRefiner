# This rule finds the extension for 3prime partial trasnripts

rule three_prime_extension:
    input:
        expand("{three_prime_extensions}/three_prime_extensions.bed",three_prime_extensions=three_prime_extensions)

rule do_three_prime_extension:
    input:
        expand("{transcript_location_in_genome}/3prime_partial_transcripts_in_genome.txt",transcript_location_in_genome=transcript_location_in_genome)
    output:
        expand("{three_prime_extensions}/three_prime_extensions.bed",three_prime_extensions=three_prime_extensions)
    threads:
        20
    shell:
        """bash bash_scripts/three_prime/3pp_RPF_extensions_finder_intron_length_200.sh {input} > {output}   """
