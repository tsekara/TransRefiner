rule BLAST_RPF_internal_extended_fasta_vs_Uniprot:
     input:
        expand("{blastx_output_dir}/Grouped_blastx_RPF_internal_extended_swissprot_trembl_max_target_3_evalue_0.001.txt",
                                                            blastx_output_dir=blastx_output_dir,
                                                            max_blast_targets=max_blast_targets,
                                                            blast_min_eval=blast_min_eval,
                                                            transcriptome_fasta_file=transcriptome_fasta_file,
                                                            Uniprot_dir=Uniprot_dir)

rule do_blastX_internal_extended_fasta:
    input:
        fasta=expand("{internal_extensions}/Extended_internal_extensions.fasta",internal_extensions=internal_extensions),
        protein_DB=expand("{Uniprot_dir}/Blastx_Uniprot_hits_seqs",
                                                            Uniprot_dir=Uniprot_dir)
    output:
        expand("{blastx_output_dir}/blastx_RPF_internal_extended_swissprot_trembl_max_target_3_evalue_0.001.txt",
                                                            blastx_output_dir=blastx_output_dir)
    message:
        "BLASTX of RPF internal extended transcripts against Uniprot Protein sequences of original transcripts"
    threads:
        20
    params:
        max_target_seqs=3,
        eval=0.001
    shell:
        """blastx -query {input.fasta} \
        -evalue {params.eval} \
        -outfmt 6 \
        -num_threads {threads} \
        -max_target_seqs {params.max_target_seqs} \
        -db {input.protein_DB} > {output} """

rule do_Group_BLASThits_internal_extended_fasta:
    input:
            blastx_internal_extended_output=expand("{blastx_output_dir}/blastx_RPF_internal_extended_swissprot_trembl_max_target_3_evalue_0.001.txt",
                                                                blastx_output_dir=blastx_output_dir),
            Uniprot_blasthits_fasta=expand("{Uniprot_dir}/uniprot_sprot_TREMBL.fasta",
                                                                Uniprot_dir=Uniprot_dir),
            internal_extended_fasta=expand("{internal_extensions}/Extended_internal_extensions.fasta",internal_extensions=internal_extensions)
    output:
            expand("{blastx_output_dir}/Grouped_blastx_RPF_internal_extended_swissprot_trembl_max_target_3_evalue_0.001.txt",
                                                                blastx_output_dir=blastx_output_dir)
    message:
            "Grouping the segmented BLASTX hits"
    shell:
            """ perl ~/Documents/trinityrnaseq-v2.10.0/util/misc/blast_outfmt6_group_segments.pl \
            {input.blastx_internal_extended_output} {input.internal_extended_fasta} {input.Uniprot_blasthits_fasta} > {output} """
