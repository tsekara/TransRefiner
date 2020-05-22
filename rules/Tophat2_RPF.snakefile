rule tophat_RPF:
        input:
            expand("{Tophat2_output}/Smed_Axl_bowtie",Tophat2_output=Tophat2_output),
            expand("{Tophat2_RPF_bam}",Tophat2_RPF_bam=Tophat2_RPF_bam)

rule do_bowtie_build:
    input:
            expand("{genome}",genome=genome)
    output:
            touch(expand("{Tophat2_output}/Smed_Axl_bowtie",Tophat2_output=Tophat2_output))
    message:
            "Bowtie2 index of Genome"
    threads:
            20
    shell:
            """bowtie-build -f {input} {output}"""


rule do_tophat_RPF:
    input:
            fastq=expand("{RPF_pooled_fastq_file}",RPF_pooled_fastq_file=RPF_pooled_fastq_file),
            index=expand("{Tophat2_output}/Smed_Axl_bowtie",Tophat2_output=Tophat2_output)
    output:
            touch(expand("{Tophat2_RPF_bam}",Tophat2_RPF_bam=Tophat2_RPF_bam))
    message:
            "Topaht2 mRNA alignment"
    threads:
            20
    params:
            g=10
    shell:
            """tophat -o {output} -p {threads} -g {params} --bowtie1 {input.index} {input.fastq} """
