rule gmap:
  input:
        expand("{Gmap_db}/{gmap_db_name}.contig",Gmap_db=Gmap_db,gmap_db_name=gmap_db_name),
        expand("{Gmap_dir}/GMAP_vs_{fasta_name}_min_idn_0.7_min_cov_0.7.gff3",
                                                            Gmap_dir=Gmap_dir,
                                                            fasta_name=fasta_name,
                                                            gmap_db_name=gmap_db_name)

rule do_gmap_build:
    input:
        expand("{genome}",genome=genome)
    output:
        touch(expand("{Gmap_db}/{gmap_db_name}.contig",Gmap_db=Gmap_db,gmap_db_name=gmap_db_name))
    shell:
        """ gmap_build -d {gmap_db_name} -D . {input} """



rule do_gmap:
    input:
        genome_fasta=expand('{genome}',genome=genome),

        transcriptome_fasta=expand('{transcriptome_dir}/{fasta_name}',
                                    transcriptome_dir=transcriptome_dir,
                                                  fasta_name=fasta_name),
        gmapdb = expand("{Gmap_db}",Gmap_db=Gmap_db),

        gmap_db_name=expand("{gmap_db_name}",gmap_db_name=gmap_db_name)

    output:
        touch(expand("{Gmap_dir}/GMAP_vs_{fasta_name}_min_idn_0.7_min_cov_0.7.gff3",
                                                            Gmap_dir=Gmap_dir,
                                                            fasta_name=fasta_name))

    params:
        min_coverage = 0.7,
        min_identity = 0.7,
        cores = 10,
        gff3format = 2

    #log:
    #    expand("{Gmap_dir}/Gmap_log.txt",Gmap_dir=Gmap_dir)

    shell:
        """
        gmap -d {input.gmap_db_name} -D {input.gmapdb} \
        -t {params.cores} -f {params.gff3format} -n1 \
        --min-trimmed-coverage {params.min_coverage} --min-identity {params.min_identity} \
        {input.transcriptome_fasta} > {output}
        """
