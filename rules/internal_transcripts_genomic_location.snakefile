# This locates the genomics location for all the incomplete internal transcripts including its
# genomic and transcriptomic strand incomplete_ids_orientation

rule internal_transcripts_genomic_location:
    input:
        expand("{transcript_location_in_genome}/internal_transcripts_in_genome.txt",transcript_location_in_genome=transcript_location_in_genome)

rule do_internal_transcripts_genomic_location:
    input:
        internal_pp=expand("{Gmap_dir}/GMAP_internal.{extn}",extn="gff3",Gmap_dir=Gmap_dir),
        internal_pp_ids_strand=expand("{Gmap_dir}/internal_ids_transcript_orientation.txt",Gmap_dir=Gmap_dir)
    output:
        internal_loc = expand("{transcript_location_in_genome}/internal_transcripts_in_genome.txt",transcript_location_in_genome=transcript_location_in_genome)
    threads:
        20
    shell:
            """
              grep "gene" {input.internal_pp} | cut -f1,4,5,7,9| awk -F"=" '$4=$1' OFS="\t" | cut -f1,2,3,4,7 | awk '{{print $5,$1,$2,$3,$4}}' OFS="\t" |  awk 'NR==FNR {{a[$1] = $2; next}} {{ $6=a[$1] ;print $0}}' OFS="\t" {input.internal_pp_ids_strand}  - | awk '{{print $2,$3,$4,$5,$1,$6}}' OFS="\t" > {output.internal_loc}
            """
