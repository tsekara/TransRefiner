rule transdecoder:
  input:
    expand("{transdecoder_dir}/{transdecoder_dir}.transdecoder_dir/{fasta_name}/{fasta_name}.transdecoder.{extn}",
            fasta_name=fasta_name,
            transcriptome_dir=transcriptome_dir,
            extn=["bed","cds","pep","gff3"],
            transdecoder_dir=transDecoder_output,
            Pfam_dir=Pfam_dir,
            Uniprot_dir=Uniprot_dir)

rule do_transdecoder:
  input:
    fasta=expand('{transcriptome_dir}/{fasta_name}',
            transcriptome_dir=transcriptome_dir,
            fasta_name=fasta_name),

    pfam=expand("{Pfam_dir}/pfam_{fasta_name}.domtblout",
            Pfam_dir=Pfam_dir,
            fasta_name=fasta_name),

    blast=expand("{Uniprot_dir}/blastx_{fasta_name}_max_target_1_evalue_1e-3.txt",
                                                        Uniprot_dir=Uniprot_dir,
                                                        fasta_name=fasta_name)
  output:
    touch(expand("{transdecoder_dir}/{fasta_name}/{fasta_name}.transdecoder.{extn}",
            fasta_name=fasta_name,
            extn=["bed","cds","pep","gff3"],
            transdecoder_dir=transDecoder_output))
  threads:
    1
  shell:
    """
    mkdir -p {transDecoder_output} && cd {transDecoder_output}
    TransDecoder.LongOrfs -t {input.fasta}
    TransDecoder.Predict -t {input.fasta} --cpu 8 --single_best_orf --retain_long_orfs 300 \
    --retain_pfam_hits {input.pfam} --retain_blastp_hits {input.blast}
    """
