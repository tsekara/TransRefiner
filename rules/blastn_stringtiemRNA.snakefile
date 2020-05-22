rule blastN_stringTie_mRNA:
        input:
            expand("{Novel_transcripts_dir}/BlastN_StringTie_mRNA_vs_dd_Smed_v6_{max_blast_targets}_evalue_{blast_min_eval}.txt",
                                                Novel_transcripts_dir=Novel_transcripts_dir,
                                                max_blast_targets=max_blast_targets,
                                                blast_min_eval=blast_min_eval)
rule do_blastN_stringTie_mRNA:
    input:
            mRNA_fasta=expand("{Novel_transcripts_dir}/mRNA_stringtie.fasta",
                                                Novel_transcripts_dir=Novel_transcripts_dir),
            dd_Smed_v6_trans_db=expand("{transcriptome_dir}/{transcriptome_db_title}",
                                                transcriptome_dir=transcriptome_dir,
                                                transcriptome_db_title=transcriptome_db_title),
                                                
    output:
            touch(expand("{Novel_transcripts_dir}/BlastN_StringTie_mRNA_vs_dd_Smed_v6_{max_blast_targets}_evalue_{blast_min_eval}.txt",
                                                            Novel_transcripts_dir=Novel_transcripts_dir,
                                                            max_blast_targets=max_blast_targets,
                                                            blast_min_eval=blast_min_eval))
    message:
            "Running BLASTN for {input.mRNA_fasta} against {input.dd_Smed_v6_trans_db}"
    threads:
            20
    params:
            no_of_targets=1
    shell:
            """blastn -query {input.mRNA_fasta} \
            -evalue 1e-3 \
            -outfmt 6 \
            -num_threads {threads} \
            -max_target_seqs {params.no_of_targets} \
            -db {input.dd_Smed_v6_trans_db} > {output} """
