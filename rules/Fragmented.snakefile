rule fragmented:
  input:
        expand("{Fragmented_transcripts_dir}/only_repeated_scaffolds.bed",
                                                            Fragmented_transcripts_dir=Fragmented_transcripts_dir)

rule do_fragmented:
    input:
        planariode_seqs=expand("{Gmap_dir}/GMAP_vs_planariode_seqs_min_idn_0.7_min_cov_0.7.gff3",
                                                            Gmap_dir=Gmap_dir),
        dd_Smed_v6_seqs=expand("{Gmap_dir}/GMAP_vs_{fasta_name}_min_idn_0.7_min_cov_0.7.gff3",
                                                            Gmap_dir=Gmap_dir,
                                                            fasta_name=fasta_name)

    output:
        touch(expand("{Fragmented_transcripts_dir}/only_repeated_scaffolds.bed",
                                                            Fragmented_transcripts_dir=Fragmented_transcripts_dir))
    shell:
        """
        bedtools intersect -a <(grep "gene" {input.dd_Smed_v6_seqs} \
         | cut -f1,4,5,9 | awk -F"=" '$1=$1' OFS="\t" - | cut -f1,2,3,6) \
          -b <(grep "gene" {input.planariode_seqs} | cut -f1,4,5 \
           | sort-bed -) | sort-bed - | perl -pe 's/_[0-9]*\w$//gm' | \
           awk '{{$5=$4; print $0}}' OFS="\t" - | \
           sort-bed -| awk '!seen[$5]++' - | cut -f1,2,3,4  \
           awk '{if (x[$1]) {x_count[$1]++; print $0; if (x_count[$1] == 1) {print x[$1]}} x[$1] = $0}' \
            -  | sort-bed> {output}
        """
