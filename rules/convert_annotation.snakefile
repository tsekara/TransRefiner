singularity: config['htseq_clip_singularity_image']

rule convert_annotation:
	input:
		expand("{annotation_dir}/{gff_name}_flattened.txt.gz", annotation_dir = annotation_dir, gff = gff_file, gff_name = gff_name)

rule do_convert_annotation:
	input:
		expand("{annotation_dir}/{gff_name}_flattened.bed.gz", annotation_dir = annotation_dir, gff = gff_file, gff_name = gff_name)
	output:
		expand("{annotation_dir}/{gff_name}_flattened.txt.gz", annotation_dir = annotation_dir, gff = gff_file, gff_name = gff_name)
	log:
		expand("{log_dir}/{gff_name}_flattened.txt.gz.log", log_dir = log_dir, gff_name = gff_name)
	message:
		"converting {input}"
	shell:
		"""mkdir -p {annotation_dir}  && mkdir -p {log_dir} && htseq-clip mapToId -a {input} -o {output} 2> {log}"""
