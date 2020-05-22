rule get_UniprotDB:
        input:
            expand("{Uniprot_dir}/uniprot_sprot_TREMBL.fasta.fai",Uniprot_dir=Uniprot_dir)

rule download_Trembl:
        output:
            expand("{Uniprot_dir}/uniprot_TREMBL.fasta.gz",Uniprot_dir=Uniprot_dir),
        message:
            "Downloading the knowledgebase-2017-10 database "
        shell:
            """ curl ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz -o {output}
            """
rule download_Swissprot:
        output:
            expand("{Uniprot_dir}/uniprot_sprot.fasta.gz",Uniprot_dir=Uniprot_dir)
        message:
            "Downloading the Swissprot-2017-10 database"
        shell:
            """ curl ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz -o {output}
            """
rule Merge_Swissprot_Trembl:
        input:
            SP=expand("{Uniprot_dir}/uniprot_sprot.fasta.gz",Uniprot_dir=Uniprot_dir),
            TR=expand("{Uniprot_dir}/uniprot_TREMBL.fasta.gz",Uniprot_dir=Uniprot_dir)
        output:
            expand("{Uniprot_dir}/uniprot_sprot_TREMBL.fasta",Uniprot_dir=Uniprot_dir)
        message:
            "Merging Swissprot and Trembl database"
        shell:
            """ cat {input.SP} {input.TR} | pigz -d > {output}
            """
rule Index_Uniprot:
        input:
            expand("{Uniprot_dir}/uniprot_sprot_TREMBL.fasta",Uniprot_dir=Uniprot_dir)
        output:
            expand("{Uniprot_dir}/uniprot_sprot_TREMBL.fasta.fai",Uniprot_dir=Uniprot_dir)
        message:
            "Indexing the fasta file"
        shell:
            """ samtools faidx {input}
            """
