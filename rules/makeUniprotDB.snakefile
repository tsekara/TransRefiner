rule makeblast_UniprotDB:
        input:
          expand("{Uniprot_dir}/Uniprot_KB_SP_seqs",  Uniprot_dir=Uniprot_dir)

rule do_makeblast_Uniprot_db:
        input:
            expand("{Uniprot_dir}/uniprot_sprot_TREMBL.fasta",Uniprot_dir=Uniprot_dir)
        output:
            touch(expand("{Uniprot_dir}/Uniprot_KB_SP_seqs",
                                Uniprot_dir=Uniprot_dir))
        message:
            "Making the BLAST database"
        params:
            Uniprot_db_title=Uniprot_db_title
        threads:
            20
        shell:
            """makeblastdb -in {input} \
             -input_type fasta \
             -dbtype prot \
             -title {params.Uniprot_db_title}  \
             -out {output} """
