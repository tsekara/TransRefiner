# This locates the genomics location for all the 3pp incomplete transcripts including its
# genomic and transcriptomic strand incomplete_ids_orientation
rule three_pp_transcripts_genomic_location:
    input:
        expand("{transcript_location_in_genome}/3prime_partial_transcripts_in_genome.txt",transcript_location_in_genome=transcript_location_in_genome),

rule do_three_pp_transcripts_genomic_location:
    input:
        three_pp=expand("{Gmap_dir}/GMAP_3prime_partial.{extn}",extn="gff3",Gmap_dir=Gmap_dir),
        three_pp_ids_strand=expand("{Gmap_dir}/3prime_partial_ids_transcript_orientation.txt",Gmap_dir=Gmap_dir),

    output:
        three_prime_loc = expand("{transcript_location_in_genome}/3prime_partial_transcripts_in_genome.txt",transcript_location_in_genome=transcript_location_in_genome),
    threads:
        20
    shell:
            """
            grep "gene" {input.three_pp} | cut -f1,4,5,7,9| awk -F"=" '$4=$1' OFS="\t" | cut -f1,2,3,4,7 | awk '{{print $5,$1,$2,$3,$4}}' OFS="\t" |  awk 'NR==FNR {{a[$1] = $2; next}} {{ $6=a[$1] ;print $0}}' OFS="\t" {input.three_pp_ids_strand}  - | awk '{{print $2,$3,$4,$5,$1,$6}}' OFS="\t" > {output.three_prime_loc}
            """
