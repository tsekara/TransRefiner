rule gmap_planariode:
  input:
        expand("{Gmap_dir}/GMAP_vs_planariode_seqs_min_idn_0.7_min_cov_0.7.gff3",
                                                            Gmap_dir=Gmap_dir,
                                                            planariode_fasta=planariode_fasta)

rule do_gmap_planariode:
    input:
        genome_fasta=expand("{genome}",genome=genome),
        planariode_fasta=expand("{planariode_fasta}",planariode_fasta=planariode_fasta),
        gmapdb = expand("{Gmap_db}",Gmap_db=Gmap_db),
        gmap_db_name=expand("{gmap_db_name}",gmap_db_name=gmap_db_name)

    output:
        touch(expand("{Gmap_dir}/GMAP_vs_planariode_seqs_min_idn_0.7_min_cov_0.7.gff3",
                                                            Gmap_dir=Gmap_dir))
    params:
        min_coverage = 0.7,
        min_identity = 0.7,
        cores = 10,
        gff3format = 2

    message:
        "Aligning the planariode seqs onto Genome"

    shell:
        """
        gmap -d {input.gmap_db_name} -D {input.gmapdb} \
        -t {params.cores} -f {params.gff3format} -n1 \
        --min-trimmed-coverage {params.min_coverage} --min-identity {params.min_identity} \
        {input.planariode_fasta} > {output}
        """
