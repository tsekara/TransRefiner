rule blastX_interpro_novel_transcripts:
        input:
            expand("{Novel_transcripts_dir}/Grouped_Blastx_Novel_transcripts_{max_blast_targets}_evalue_{blast_min_eval}.txt",
                                                            Novel_transcripts_dir=Novel_transcripts_dir,
                                                            max_blast_targets=max_blast_targets,
                                                            blast_min_eval=blast_min_eval),
            expand("{Novel_transcripts_dir}/Final_novel_transcripts.gff3",
                                                            Novel_transcripts_dir=Novel_transcripts_dir

rule do_blastX_novel_transcripts:
    input:
            fasta=expand("{Novel_transcripts_dir}/Final_novel_transcripts.fasta",
                                                            Novel_transcripts_dir=Novel_transcripts_dir),
            protein_DB=expand("{Uniprot_dir}/{Uniprot_db_title}",
                                                            Uniprot_dir=Uniprot_dir,
                                                            Uniprot_db_title=Uniprot_db_title)
    output:
            touch(expand("{Novel_transcripts_dir}/Blastx_Novel transcripts_max_target_{max_blast_targets}_evalue_{blast_min_eval}.txt",
                                                            blastx_output_dir=blastx_output_dir,
                                                            max_blast_targets=max_blast_targets,
                                                            blast_min_eval=blast_min_eval))
    message:
            "Running BLASTX for {input.fasta} against {input.protein_DB}"
    threads:
            20
    params:
            no_of_targets=1
    shell:
            """blastx -query {input.fasta} \
            -evalue 1e-3 \
            -outfmt 6 \
            -num_threads {threads} \
            -max_target_seqs {params.no_of_targets} \
            -db {input.protein_DB} > {output} """

rule do_group_blastHits_novel_transcripts:
    input:
            blastx_output=expand("{Novel_transcripts_dir}/Blastx_Novel transcripts_max_target_{max_blast_targets}_evalue_{blast_min_eval}.txt",
                                                            blastx_output_dir=blastx_output_dir,
                                                            max_blast_targets=max_blast_targets,
                                                            blast_min_eval=blast_min_eval),
            Uniprot_fasta=expand("{Uniprot_dir}/uniprot_sprot_TREMBL.fasta",Uniprot_dir=Uniprot_dir),
            fasta=expand("{Novel_transcripts_dir}/Final_novel_transcripts.fasta",
                                                            Novel_transcripts_dir=Novel_transcripts_dir)
    output:
            touch(expand("{Novel_transcripts_dir}/Grouped_Blastx_Novel_transcripts_{max_blast_targets}_evalue_{blast_min_eval}.txt",
                                                            Novel_transcripts_dir=Novel_transcripts_dir,
                                                            max_blast_targets=max_blast_targets,
                                                            blast_min_eval=blast_min_eval))
    message:
            "Grouping the segmented BLASTX hits"
    shell:
            """perl ~/Documents/trinityrnaseq-v2.10.0/util/misc/blast_outfmt6_group_segments.pl \
            {input.blastx_output} {input.fasta} {input.Uniprot_fasta} > {output}"""

rule do_interpro:
    input:
            fasta=expand("{Novel_transcripts_dir}/Final_novel_transcripts.fasta",
                                                            Novel_transcripts_dir=Novel_transcripts_dir)
    output:
            touch(expand("{Novel_transcripts_dir}/Final_novel_transcripts.gff3",
                                                            Novel_transcripts_dir=Novel_transcripts_dir))
    message:
            "Grouping the segmented BLASTX hits"
    shell:
            """bash ~/ownCloud/conda_xtender/interproscan.sh  \
                -cpu 10 \
                -f tsv,gff3 \
              -goterms -i {input} \
              -pa -t n -d {output} """
