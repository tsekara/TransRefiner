rule pfam:
        input:
            expand("{Pfam_dir}/pfam_{fasta_name}.domtblout",
                    Pfam_dir=Pfam_dir,
                    fasta_name=fasta_name)

rule do_pfam:
        input:
            fasta=expand('{transcriptome_dir}/{fasta_name}',
                            transcriptome_dir=transcriptome_dir,
                            fasta_name=fasta_name),

            hmm_db=expand("{Pfam_hmm}",Pfam_hmm=Pfam_hmm)
        output:
            expand("{Pfam_dir}/pfam_{fasta_name}.domtblout",
                    Pfam_dir=Pfam_dir,
                    fasta_name=fasta_name)
        threads:
            20
        shell:
            """hmmscan --cpu {threads} \
            --domtblout {output} \
            {input.hmm_db} \
            {input.fasta} """
