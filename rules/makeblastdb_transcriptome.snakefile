rule makeblast_dd_Smed_v6_transcriptome:
        input:
          expand("{transcriptome_dir}/{transcriptome_db_title}",
                              transcriptome_dir=transcriptome_dir,
                              transcriptome_db_title=transcriptome_db_title)

rule do_dd_Smed_v6_transcriptome:
        input:
            expand("{transcriptome_dir}/{transcriptome_fasta}",
                                        transcriptome_dir=transcriptome_dir,
                                        transcriptome_fasta=transcriptome_fasta)
        output:
            touch(expand("{transcriptome_dir}/{transcriptome_db_title}",
                                transcriptome_dir=transcriptome_dir,
                                transcriptome_db_title=transcriptome_db_title))
        message:
            "Making the BLAST database for {transcriptome_fasta}"
        params:
            transcriptome_db_title=transcriptome_db_title
        threads:
            20
        shell:
            """makeblastdb -in {input} \
             -input_type fasta \
             -dbtype nucl \
             -title {params.transcriptome_db_title}  \
             -out {output} """
