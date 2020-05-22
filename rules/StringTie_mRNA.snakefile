rule mRNA_stringtie:
        input:
          expand("{Novel_transcripts_dir}/mRNA_stringtie.fasta",Novel_transcripts_dir=Novel_transcripts_dir)

rule do_mRNA_stringtie:
        input:
            mRNA_tophat_bam=expand("{Tophat2_mRNA_bam}",Tophat2_mRNA_bam=Tophat2_mRNA_bam),
            gmap_annotation=expand("{Gmap_dir}/GMAP_vs_{fasta_name}_min_idn_0.7_min_cov_0.7.gff3",
                                                                Gmap_dir=Gmap_dir,
                                                                fasta_name=fasta_name)
        output:
            touch(expand("{Novel_transcripts_dir}/mRNA_stringtie.gtf",Novel_transcripts_dir=Novel_transcripts_dir))
        message:
            "mRNA: stringTie assembly"
        params:
            p=20,
            C="/Users/sekaran/ownCloud/conda_xtender/Novel_transcripts//mRNA_known_trans",
            f=0.5
        threads:
            20
        shell:
            """ stringtie {input.mRNA_tophat_bam} \
             -G {input.gmap_annotation} \
             -p {params.p} \
             -o {output} \
             -v \
             -f {params.f} \
             -C {params.C} """

rule do_mRNA_convert_gtf_to_fa:
        input:
            mRNA_gtf=expand("{Novel_transcripts_dir}/mRNA_stringtie.gtf",Novel_transcripts_dir=Novel_transcripts_dir),
            genome=expand("{genome}",genome=genome)
        output:
            touch(expand("{Novel_transcripts_dir}/mRNA_stringtie.fasta",Novel_transcripts_dir=Novel_transcripts_dir))
        message:
            "mRNA: Converting the mRNA GTF into fasta file"
        threads:
            20
        shell:
            """ gffread {input.mRNA_gtf} -g {input.genome} -w {output}"""
