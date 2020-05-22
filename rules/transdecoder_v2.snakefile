rule all:
  input:
    "test.fasta.transdecoder_dir/longest_orfs.pep",
    "test.fasta.transdecoder.pep"

rule transdecoder_longORFs_predict:
  input:
    "test.fasta"
  output:
    "test.fasta.transdecoder_dir/longest_orfs.pep"
  threads:
    20
  shell:
    """
    TransDecoder.LongOrfs -t {input}
    """

rule transdecoder_predict:
    input:
        'test.fasta'
    output:
        'test.fasta.transdecoder.pep'
    threads:
        20
    shell:
        """ TransDecoder.Predict -t {input} \
        --cpu 30 --single_best_orf \
        --retain_long_orfs 300 \
        --retain_pfam_hits Pfam/pfam_dd_Smed_v6.pcf.contigs.domtblout \
        --retain_blastp_hits Uniprot_db/Blastx_dd_Smed_v6_pcf_contigs_vs_Swissprot_Trembl_Oct2017_max_target_1_evalue_1e-3.txt """
