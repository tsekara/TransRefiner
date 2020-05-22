rule incomplete_transcripts_location:
    input:
        expand("{Gmap_dir}/GMAP_5prime_partial.{extn}",extn="gff3",Gmap_dir=Gmap_dir),
        expand("{Gmap_dir}/GMAP_3prime_partial.{extn}",extn="gff3",Gmap_dir=Gmap_dir),
        expand("{Gmap_dir}/GMAP_internal.{extn}",extn="gff3",Gmap_dir=Gmap_dir)

rule do_incomplete_transcripts_location:
    input:
        GMAP_file=expand("{Gmap_dir}/GMAP_vs_{fasta_name}_min_idn_0.7_min_cov_0.7.gff3",
                                                        Gmap_dir=Gmap_dir,
                                                        fasta_name=fasta_name,
                                                        gmap_db_name=gmap_db_name),

        transdecoder_file=expand("{transdecoder_dir}/{fasta_name}/{fasta_name}.transdecoder.{extn}",
                                                extn="gff3",
                                                transdecoder_dir=transDecoder_output,
                                                fasta_name=fasta_name)
    output:
        five_pp=expand("{Gmap_dir}/GMAP_5prime_partial.{extn}",extn="gff3",Gmap_dir=Gmap_dir),
        three_pp=expand("{Gmap_dir}/GMAP_3prime_partial.{extn}",extn="gff3",Gmap_dir=Gmap_dir),
        internal=expand("{Gmap_dir}/GMAP_internal.{extn}",extn="gff3",Gmap_dir=Gmap_dir),

    threads:
        20

    #shell:
    #    """grep -w -F -f <(grep  "gene" {input.transdecoder_file}| grep "5prime_partial" - | cut -f1) {input.GMAP_file} > {output.five_pp}
    #    grep -w -F -f <(grep  "gene" {input.transdecoder_file}|grep "3prime_partial" - | cut -f1) {input.GMAP_file} > {output.three_pp}
    #    grep -w -F -f <(grep  "gene" {input.transdecoder_file}|grep "internal" - | cut -f1) {input.GMAP_file} > {output.internal}"""
    run:
        shell("""grep -w -F -f <(grep  "gene" {input.transdecoder_file}| grep "5prime_partial" - | cut -f1) {input.GMAP_file} > {output.five_pp}""")
        shell("""grep -w -F -f <(grep  "gene" {input.transdecoder_file}|grep "3prime_partial" - | cut -f1) {input.GMAP_file} > {output.three_pp}""")
        shell("""grep -w -F -f <(grep  "gene" {input.transdecoder_file}|grep "internal" - | cut -f1) {input.GMAP_file} > {output.internal}""")
