rule final_novel_transcripts:
        input:
            expand("{Novel_transcripts_dir}/Final_novel_transcripts.fasta",
                                                            Novel_transcripts_dir=Novel_transcripts_dir)
rule do_novel_trans_mRNA:
    input:
            mRNA_stringtie_fa=expand("{Novel_transcripts_dir}/mRNA_stringtie.fasta",
                                                            Novel_transcripts_dir=Novel_transcripts_dir),
            BlastN_mRNA_stringtie=expand("{Novel_transcripts_dir}/BlastN_StringTie_mRNA_vs_dd_Smed_v6_{max_blast_targets}_evalue_{blast_min_eval}.txt",
                                                            Novel_transcripts_dir=Novel_transcripts_dir,
                                                            max_blast_targets=max_blast_targets,
                                                            blast_min_eval=blast_min_eval),
            mRNA_stringtie_gtf=expand("{Novel_transcripts_dir}/mRNA_stringtie.gtf",
                                                            Novel_transcripts_dir=Novel_transcripts_dir)
    output:
            touch(expand("{Novel_transcripts_dir}/mRNA_stringtie_non_dd_Smed_v6_trasncripts.gtf",
                                                            Novel_transcripts_dir=Novel_transcripts_dir))
    message:
            "Retaining only the mRNA-based non-dd_Smed_v6 transcripts"
    threads:
            20
    shell:
            """ grep -w -F -f <(comm -23 <(grep ">" {input.mRNA_stringtie_fa} | perl -p -e 's/>//g' | sort) \
               <(cut -f1 {input.BlastN_mRNA_stringtie} | sort) | uniq) {input.mRNA_stringtie_gtf} \
               | grep -v -w "reference_id" > {output} """


rule do_novel_trans_RPF:
    input:
            RPF_stringtie_fa=expand("{Novel_transcripts_dir}/RPF_stringtie.fasta",
                                                            Novel_transcripts_dir=Novel_transcripts_dir),
            BlastN_RPF_stringtie=expand("{Novel_transcripts_dir}/BlastN_StringTie_RPF_vs_dd_Smed_v6_{max_blast_targets}_evalue_{blast_min_eval}.txt",
                                                            Novel_transcripts_dir=Novel_transcripts_dir,
                                                            max_blast_targets=max_blast_targets,
                                                            blast_min_eval=blast_min_eval),
            RPF_stringtie_gtf=expand("{Novel_transcripts_dir}/mRNA_stringtie.gtf",
                                                            Novel_transcripts_dir=Novel_transcripts_dir)
    output:
            touch(expand("{Novel_transcripts_dir}/RPF_stringtie_non_dd_Smed_v6_trasncripts.gtf", Novel_transcripts_dir=Novel_transcripts_dir))
    message:
            "Retaining only the RPF-based non-dd_Smed_v6 transcripts"
    threads:
            20
    shell:
            """ grep -w -F -f <(comm -23 <(grep ">" {input.RPF_stringtie_fa} | perl -p -e 's/>//g' | sort) \
               <(cut -f1 {input.BlastN_RPF_stringtie} | sort) | uniq) {input.RPF_stringtie_gtf} \
               | grep -v -w "reference_id" > {output} """

rule do_intersect_mRNA_RPF:
    input:
            RPF_uniq=expand("{Novel_transcripts_dir}/RPF_stringtie_non_dd_Smed_v6_trasncripts.gtf",
                                                            Novel_transcripts_dir=Novel_transcripts_dir),
            mRNA_uniq=expand("{Novel_transcripts_dir}/mRNA_stringtie_non_dd_Smed_v6_trasncripts.gtf",
                                                            Novel_transcripts_dir=Novel_transcripts_dir)
    output:
            touch(expand("{Novel_transcripts_dir}/Final_novel_trasncripts_mRNA_RPF_coverage.bed",
                                                            Novel_transcripts_dir=Novel_transcripts_dir))
    message:
            "Retaining only the overalpping mRNA-based and RPF-based non-dd_Smed_v6 transcripts"
    threads:
            20
    shell:
            """ bedtools intersect -wa -wb -a <( grep -w "transcript" {input.mRNA_uniq} | cut -f1,4,5,9) \
             -b <( grep -w "transcript" {input.RPF_uniq} | cut -f1,4,5,9) \
              > {output} """

rule do_extract_novel_mRNAs_with_RPF_coverage:
    input:
            mRNA_RPF_overlap=expand("{Novel_transcripts_dir}/Final_novel_trasncripts_mRNA_RPF_coverage.bed",
                                                            Novel_transcripts_dir=Novel_transcripts_dir),
            mRNA_stringtie_gtf=expand("{Novel_transcripts_dir}/mRNA_stringtie.gtf",
                                                            Novel_transcripts_dir=Novel_transcripts_dir)
    output:
            touch(expand("{Novel_transcripts_dir}/Final_novel_trasncripts_mRNA_stringtie.gtf",
                                                            Novel_transcripts_dir=Novel_transcripts_dir))
    message:
            "extracting the transcript-ids of overlapping intervals and extracting the mRNA-stringtie gtf entries"
    threads:
            20
    shell:
            """ grep -w -F -f <(awk '$1=$1' FS=";" OFS="\t" {input.mRNA_RPF_overlap} \
             | cut -f5 | perl -p -e 's/ transcript_id |"//g') {input.mRNA_stringtie_gtf} \
             > {output} """

rule do_convert_gtf_to_fasta:
    input:
            genome_fasta=expand("{genome}", genome=genome),
            final_mRNA_stringtie_gtf=expand("{Novel_transcripts_dir}/Final_novel_trasncripts_mRNA_stringtie.gtf",
                                                            Novel_transcripts_dir=Novel_transcripts_dir)
    output:
            touch(expand("{Novel_transcripts_dir}/Final_novel_transcripts.fasta",
                                                            Novel_transcripts_dir=Novel_transcripts_dir))
    message:
            "extracting the transcript-ids of overlapping intervals and extracting the mRNA-stringtie gtf entries"
    threads:
            20
    shell:
            """ gffread {input.final_mRNA_stringtie_gtf} \
             -g {input.genome_fasta} \
             -w {output} """
