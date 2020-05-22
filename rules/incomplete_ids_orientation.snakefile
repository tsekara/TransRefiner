rule incomplete_transcripts_ids_orientation:
    input:
        expand("{Gmap_dir}/5prime_partial_ids_transcript_orientation.txt",Gmap_dir=Gmap_dir),
        expand("{Gmap_dir}/3prime_partial_ids_transcript_orientation.txt",Gmap_dir=Gmap_dir),
        expand("{Gmap_dir}/internal_ids_transcript_orientation.txt",Gmap_dir=Gmap_dir),

rule do_incomplete_transcripts_ids_orientation:
    input:
        five_pp=expand("{Gmap_dir}/GMAP_5prime_partial.{extn}",extn="gff3",Gmap_dir=Gmap_dir),
        three_pp=expand("{Gmap_dir}/GMAP_3prime_partial.{extn}",extn="gff3",Gmap_dir=Gmap_dir),
        internal=expand("{Gmap_dir}/GMAP_internal.{extn}",extn="gff3",Gmap_dir=Gmap_dir),

        transdecoder_file=expand("{transdecoder_dir}/{fasta_name}/{fasta_name}.transdecoder.{extn}",
                                                extn="gff3",
                                                transdecoder_dir=transDecoder_output,
                                                fasta_name=fasta_name)
    output:
        five_pp=expand("{Gmap_dir}/5prime_partial_ids_transcript_orientation.txt",Gmap_dir=Gmap_dir),
        three_pp=expand("{Gmap_dir}/3prime_partial_ids_transcript_orientation.txt",Gmap_dir=Gmap_dir),
        internal=expand("{Gmap_dir}/internal_ids_transcript_orientation.txt",Gmap_dir=Gmap_dir),

    threads:
        20

    shell:
            """grep -w -F -f <(grep "gene" {input.five_pp} | awk -F'[=.]' '$1=$1' OFS="\\t" - | awk '{{print $NF}}' - OFS="\\t") <(grep "gene" {input.transdecoder_file} | awk '{{print $1"\\t"$7}}' -) > {output.five_pp}
            grep -w -F -f <(grep "gene" {input.three_pp} | awk -F'[=.]' '$1=$1' OFS="\\t" - | awk '{{print $NF}}' - OFS="\\t") <(grep "gene" {input.transdecoder_file} | awk '{{print $1"\\t"$7}}' -) > {output.three_pp}
            grep -w -F -f <(grep "gene" {input.internal} | awk -F'[=.]' '$1=$1' OFS="\\t" - | awk '{{print $NF}}' - OFS="\\t") <(grep "gene" {input.transdecoder_file} | awk '{{print $1"\\t"$7}}' -) > {output.internal}"""

#    run:
#        shell("""grep -w -F -f <(grep "gene" {input.five_pp} | awk -F'[=.]' '$1=$1' OFS="\t" - | awk '{{print $9}}' - OFS="\t") <(grep "gene" {input.transdecoder_file} | awk '{{print $1,$7}}' -  OFS="\t") > {output.five_pp}""")
#        shell("""grep -w -F -f <(grep "gene" {input.three_pp} | awk -F'[=.]' '$1=$1' OFS="\t" - | awk '{{print $9}}' - OFS="\t") <(grep "gene" {input.transdecoder_file} | awk '{{print $1,$7}}' -  OFS="\t") > {output.three_pp}""")
#        shell("""grep -w -F -f <(grep "gene" {input.internal} | awk -F'[=.]' '$1=$1' OFS="\t" - | awk '{{print $9}}' - OFS="\t") <(grep "gene" {input.transdecoder_file} | awk '{{print $1,$7}}' -  OFS="\t") > {output.internal}""")
