rule makeblast_BLAST_hits:
     input:
        expand("{Uniprot_dir}/Blastx_Uniprot_hits_seqs",
                                            Uniprot_dir=Uniprot_dir,
                                            blastx_output_dir=blastx_output_dir,
                                            max_blast_targets=max_blast_targets,
                                            blast_min_eval=blast_min_eval,
                                            transcriptome_fasta_file=transcriptome_fasta_file)

rule do_extract_UniProt_seqs_of_BLAST_hits:
     input:
         blast_output=expand("{blastx_output_dir}/Blastx_dd_Smed_v6_pcf_contigs_vs_Swissprot_Trembl_Oct2017_max_target_{max_blast_targets}_evalue_{blast_min_eval}.txt",
                                                         blastx_output_dir=blastx_output_dir,
                                                         max_blast_targets=max_blast_targets,
                                                         blast_min_eval=blast_min_eval),
                                                         Uniprot_SP_TR=expand("{Uniprot_dir}/uniprot_sprot_TREMBL.fasta", Uniprot_dir=Uniprot_dir)
     output:
         expand("{Uniprot_dir}/Blastx_hits_uniprot.fasta", Uniprot_dir=Uniprot_dir)
     threads:
         40
     message:
        "Extracting the Uniprot seqs of only BLAST hits for dd_Smed_v6 transcriptome"
     shell:
            """ cut -f2 {input.blast_output} | xargs samtools faidx {input.Uniprot_SP_TR} > {output}  """

rule do_makeblast_BLAST_hits:
        input:
            expand("{Uniprot_dir}/Blastx_hits_uniprot.fasta", Uniprot_dir=Uniprot_dir)
        output:
            touch(expand("{Uniprot_dir}/Blastx_Uniprot_hits_seqs",
                                Uniprot_dir=Uniprot_dir))
        message:
            "Making the BLASTX database for only Prot seqs of BlastX hits"
        params:
            Uniprot_db_title="Blastx_Uniprot_hits_seqs"
        threads:
            20
        shell:
            """makeblastdb -in {input} \
             -input_type fasta \
             -dbtype prot \
             -title {params.Uniprot_db_title}  \
             -out {output}"""
