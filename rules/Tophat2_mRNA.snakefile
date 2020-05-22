rule tophat_mRNA:
        input:
            expand("{Tophat2_output}/Smed_Axl_bowtie2",Tophat2_output=Tophat2_output),
            expand("{Tophat2_mRNA_bam}",Tophat2_mRNA_bam=Tophat2_mRNA_bam)

rule do_bowtie2_buid:
    input:
            expand("{genome}",genome=genome)
    output:
            touch(expand("{Tophat2_output}/Smed_Axl_bowtie2",Tophat2_output=Tophat2_output))
    message:
            "Bowtie2 index of Genome"
    threads:
            20
    shell:
            """bowtie2-build -f {input} {output}"""
rule do_tophat_mRNA:
    input:
            fastq=expand("{mRNA_pooled_fastq_file}",mRNA_pooled_fastq_file=mRNA_pooled_fastq_file),
            index=expand("{Tophat2_output}/Smed_Axl_bowtie2",Tophat2_output=Tophat2_output)
    output:
            touch(expand("{Tophat2_mRNA_bam}",Tophat2_mRNA_bam=Tophat2_mRNA_bam))
    message:
            "Topaht2 mRNA alignment"
    threads:
            20
    shell:
            """tophat -o {output} -p {threads} {input.fastq} {input.index} """
