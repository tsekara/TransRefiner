rule blastX:
        input:
            expand("{blastx_output_dir}/Grouped_Blastx_dd_Smed_v6_pcf_contigs_vs_Swissprot_Trembl_Oct2017_max_target_{max_blast_targets}_evalue_{blast_min_eval}.txt",
                                                            blastx_output_dir=blastx_output_dir,
                                                            max_blast_targets=max_blast_targets,
                                                            blast_min_eval=blast_min_eval,
                                                            transcriptome_fasta_file=transcriptome_fasta_file)
rule do_blastX:
    input:
            fasta=expand('{transcriptome_dir}/{fasta_name}',transcriptome_dir=transcriptome_dir,fasta_name=fasta_name),

            protein_DB=expand("{Uniprot_dir}/{Uniprot_db_title}",Uniprot_dir=Uniprot_dir,Uniprot_db_title=Uniprot_db_title)
    output:
            expand("{blastx_output_dir}/Blastx_dd_Smed_v6_pcf_contigs_vs_Swissprot_Trembl_Oct2017_max_target_{max_blast_targets}_evalue_{blast_min_eval}.txt",
                                                            blastx_output_dir=blastx_output_dir,
                                                            max_blast_targets=max_blast_targets,
                                                            blast_min_eval=blast_min_eval)
    message:
            "Running BLASTX for {input.fasta} against {input.protein_DB}"
    threads:
            20
    params:
            no_of_targets=10
    shell:
            """blastx -query {input.fasta} \
            -evalue 1e-3 \
            -outfmt 6 \
            -num_threads {threads} \
            -max_target_seqs {params.no_of_targets} \
            -db {input.protein_DB} > {output} """

rule do_group_blastHits:
    input:
            blastx_output=expand("{blastx_output_dir}/Blastx_dd_Smed_v6_pcf_contigs_vs_Swissprot_Trembl_Oct2017_max_target_{max_blast_targets}_evalue_{blast_min_eval}.txt",
                                                            blastx_output_dir=blastx_output_dir,
                                                            max_blast_targets=max_blast_targets,
                                                            blast_min_eval=blast_min_eval),
            Uniprot_fasta=expand("{Uniprot_dir}/uniprot_sprot_TREMBL.fasta",Uniprot_dir=Uniprot_dir)
    output:
            expand("{blastx_output_dir}/Grouped_Blastx_dd_Smed_v6_pcf_contigs_vs_Swissprot_Trembl_Oct2017_max_target_{max_blast_targets}_evalue_{blast_min_eval}.txt",
                                                            blastx_output_dir=blastx_output_dir,
                                                            max_blast_targets=max_blast_targets,
                                                            blast_min_eval=blast_min_eval)
    params:
            transcriptome_fasta_file=transcriptome_fasta_file
    message:
            "Grouping the segmented BLASTX hits"
    shell:
            """ perl ~/Documents/trinityrnaseq-v2.10.0/util/misc/blast_outfmt6_group_segments.pl \
            {input.blastx_output} {params.transcriptome_fasta_file} {input.Uniprot_fasta} > {output} """
