########################
#Xtender
#author: TS
#sekaran@embl.de

########################
description = """
Workflow for XTENDER
project dir: /g/hentze/users/Thilli/Final_PhD_analysis/conda_pipeline/
contact: Thileepan Sekaran
"""

print (description)
import glob
import os
import argparse
import yaml

DEFAULT_CONFIG_FILE = "./xtender.config.yaml"

# ----------------------------------------------------------------------------------------
# Open the config file
# ----------------------------------------------------------------------------------------

with open(DEFAULT_CONFIG_FILE,"r") as config_file:
	try:
		config =yaml.safe_load(config_file)
	except yaml.YAMLError as exc:
		print (exc)


transcriptome		= 	config["transcriptome_fasta_file"]
transcriptome_fasta_file = config["transcriptome_fasta_file"]
transcriptome_fasta	=	config["transcriptome_fasta"]
transcriptome_db_title = config["transcriptome_db_title"]
transDecoder_output	=	config["transdecoder_output_dir"]
transcriptome_dir	=	config["transcriptome_dir"]
Uniprot_dir			= 	config["Uniprot_dir"]
Uniprot_fasta		=	config["Uniprot_fasta"]
Pfam_hmm			=	config["Pfam_hmm"]
Pfam_dir			=	config["Pfam_dir"]
genome_dir			=	config["genome_dir"]
genome				=	config["genome"]
gmap_db_name		=	config["gmap_db_name"]
Gmap_db				= 	config["Gmap_db"]
Gmap_dir			= 	config["Gmap_dir"]
search_length		=	config["search_length"]
RPF_pooled_fastq_file	=	config["RPF_pooled_fastq_file"]
mRNA_pooled_fastq_file	=	config["mRNA_pooled_fastq_file"]
Tophat2_output		=	config["Tophat2_output"]
Tophat2_RPF_bam		=	config["Tophat2_RPF_bam"]
Tophat2_mRNA_bam	=	config["Tophat2_mRNA_bam"]
transcript_location_in_genome = config["transcript_location_in_genome"]
three_prime_extensions = config["three_prime_extensions"]
five_prime_extensions = config["five_prime_extensions"]
internal_extensions = config["internal_extensions"]
Uniprot_db_title = config["Uniprot_db_title"]
max_blast_targets =config["max_blast_targets"]
blast_min_eval	= config["blast_min_eval"]
blastx_output_dir = config["blastx_output_dir"]
Novel_transcripts_dir = config["Novel_transcripts_dir"]
planariode_fasta=config["planariode_fasta"]
Fragmented_transcripts_dir=config["Fragmented_transcripts_dir"]

fasta_name=os.path.basename(transcriptome)

# ----------------------------------------------------------------------------------------
# COLOURS
# ----------------------------------------------------------------------------------------
GREEN   = "\x1b[32;01m" ; CLEAN   = "\x1b[39;49;00m"; RED= "\x1b[31;01m" ; MAGENTA = "\x1b[35;01m"

# Aliging the pooled mRNA fastq files with Tophat2
# ---------------------------------------------------------------------------------------
include: os.path.join(".","rules","Tophat2_mRNA.snakefile")

# Aliging the pooled RPF fastq files with Tophat2
# ---------------------------------------------------------------------------------------
include: os.path.join(".","rules","Tophat2_RPF.snakefile")

# Downlaoding the Uniprot DB via curl or (wget)
# ---------------------------------------------------------------------------------------
include: os.path.join(".","rules","getUniprotDB.snakefile")

# makeblast Uniprot_DB
# ---------------------------------------------------------------------------------------
include: os.path.join(".","rules","makeUniprotDB.snakefile")

# Transdecoder
# ---------------------------------------------------------------------------------------
include: os.path.join(".","rules","Transdecoder.snakefile")

# BLAST transcriptome vs Uniprot SP-TR
# ---------------------------------------------------------------------------------------
include: os.path.join(".","rules","blastx_uniprot.snakefile")


# extract and makeblastdb only for blasthits
# ---------------------------------------------------------------------------------------
include: os.path.join(".","rules","makeblastdb_Uniprot_BLASThits.snakefile")


# Pfam
# ---------------------------------------------------------------------------------------
include: os.path.join(".","rules","Pfam.snakefile")

# GMAP the dd_Smed_v6
# ---------------------------------------------------------------------------------------
include: os.path.join(".","rules","gmap.snakefile")


# GMAP the planariode sequences
# ---------------------------------------------------------------------------------------
include: os.path.join(".","rules","gmap_planariode.snakefile")


# incomplete_transcripts_location
# ---------------------------------------------------------------------------------------
include: os.path.join(".","rules","GMAP_incomplete_location.snakefile")


# Extracting the incomplete transcripts and its orientation in transcriptome
# ---------------------------------------------------------------------------------------
include: os.path.join(".","rules","incomplete_ids_orientation.snakefile")


# Locating the 5pp incomplete transcripts in the genome
# ---------------------------------------------------------------------------------------
include: os.path.join(".","rules","five_prime_transcripts_genomic_location.snakefile")


# Locating the 3pp incomplete transcripts in the genome
# ---------------------------------------------------------------------------------------
include: os.path.join(".","rules","three_prime_transcripts_genomic_location.snakefile")


# Locating the internal incomplete transcripts in the genome
# ---------------------------------------------------------------------------------------
include: os.path.join(".","rules","internal_transcripts_genomic_location.snakefile")


# Finding three prime extensions
# ---------------------------------------------------------------------------------------
include: os.path.join(".","rules","three_prime_extensions.snakefile")


# Cleaning the three prime extensions
# ---------------------------------------------------------------------------------------
include: os.path.join(".","rules","Clean_3pp_extensions.snakefile")


# Cleaning the three prime extensions
# ---------------------------------------------------------------------------------------
include: os.path.join(".","rules","Append_3pp_extension.snakefile")


# Finding five prime extensions
# ---------------------------------------------------------------------------------------
include: os.path.join(".","rules","five_prime_extension.snakefile")


# Cleaning five prime extensions
# ---------------------------------------------------------------------------------------
include: os.path.join(".","rules","Clean_5pp_extensions.snakefile")


# Appending the five prime extensions to original fasta
# ---------------------------------------------------------------------------------------
include: os.path.join(".","rules","Append_5pp_extension.snakefile")



# 5pp RPF extended fasta vs Unirpot
# ---------------------------------------------------------------------------------------
include: os.path.join(".","rules","blastx_RPF_5pp_extended_Vs_Uniprot.snakefile")


# Finding the extesnions of internal transcripts
# ---------------------------------------------------------------------------------------
include: os.path.join(".","rules","internal_extensions.snakefile")


# Cleaning internal extensions
# ---------------------------------------------------------------------------------------
include: os.path.join(".","rules","Clean_internal_extensions.snakefile")


# Appending internal extensions
# ---------------------------------------------------------------------------------------
include: os.path.join(".","rules","Append_internal_extension.snakefile")


# internal RPF extended fasta vs Unirpot
# ---------------------------------------------------------------------------------------
include: os.path.join(".","rules","blastx_RPF_internal_extended_Vs_Uniprot.snakefile")



# Finding the 3pp extensions
# ---------------------------------------------------------------------------------------
include: os.path.join(".","rules","test_3pp_extensions_finder.snakefile")


# Filtering the 3pp extensions
# ---------------------------------------------------------------------------------------
include: os.path.join(".","rules","test_3pp_extensions_filter.snakefile")


# Appending the 3pp extensions
# ---------------------------------------------------------------------------------------
include: os.path.join(".","rules","test_3pp_extensions_appender.snakefile")

# 3pp RPF extended fasta vs Unirpot
# ---------------------------------------------------------------------------------------
include: os.path.join(".","rules","blastx_RPF_3pp_extended_Vs_Uniprot.snakefile")


# Finding the mRNA 3pp extensions
# ---------------------------------------------------------------------------------------
include: os.path.join(".","rules","mRNA_3primepartial_extensions_finder.snakefile")

# Filtering the mRNA 3pp extensions
# ---------------------------------------------------------------------------------------
include: os.path.join(".","rules","mRNA_3primepartial_extensions_filter.snakefile")

# Appending the mRNA 3pp extensions
# ---------------------------------------------------------------------------------------
include: os.path.join(".","rules","mRNA_3pp_extensions_appender.snakefile")


# Finding the mRNA 5pp extensions
# ---------------------------------------------------------------------------------------
include: os.path.join(".","rules","mRNA_5primepartial_extensions_finder.snakefile")


# Filtering the mRNA 5pp extensions
# ---------------------------------------------------------------------------------------
include: os.path.join(".","rules","mRNA_5primepartial_extensions_filter.snakefile")


# Appending the mRNA 5pp extensions
# ---------------------------------------------------------------------------------------
include: os.path.join(".","rules","mRNA_5primepartial_extensions_append.snakefile")


# Finding the mRNA internal extensions
# ---------------------------------------------------------------------------------------
include: os.path.join(".","rules","mRNA_internal_extensions_finder.snakefile")

# Filtering the mRNA internal extensions
# ---------------------------------------------------------------------------------------
include: os.path.join(".","rules","mRNA_internal_extensions_filter.snakefile")

# Appending the mRNA internal extensions
# ---------------------------------------------------------------------------------------
include: os.path.join(".","rules","mRNA_internal_extensions_append.snakefile")


# BLAST mRNA 5pp extended fasta vs Uniprot
# ---------------------------------------------------------------------------------------
include: os.path.join(".","rules","blastx_mRNA_5pp_extended_Vs_Uniprot.snakefile")

# BLAST mRNA 5pp extended fasta vs Uniprot
# ---------------------------------------------------------------------------------------
include: os.path.join(".","rules","blastx_mRNA_3pp_extended_Vs_Uniprot.snakefile")


# BLAST mRNA internal extended fasta vs Uniprot
# ---------------------------------------------------------------------------------------
include: os.path.join(".","rules","blastx_mRNA_internal_extended_Vs_Uniprot.snakefile")

# makeBLAST database for transcriptome
# ---------------------------------------------------------------------------------------
include: os.path.join(".","rules","makeblastdb_transcriptome.snakefile")


# stringTie for Tophat2-mRNA
# ---------------------------------------------------------------------------------------
include: os.path.join(".","rules","StringTie_mRNA.snakefile")

# BLASTN stringtie mRNA.fasta vs dd_Smed_v6.fasta
# ---------------------------------------------------------------------------------------
include: os.path.join(".","rules","blastn_stringtiemRNA.snakefile")


# stringTie for Tophat2-RPF
# ---------------------------------------------------------------------------------------
include: os.path.join(".","rules","StringTie_RPF.snakefile")


# BLASTN stringtie RPF.fasta vs dd_Smed_v6.fasta
# ---------------------------------------------------------------------------------------
include: os.path.join(".","rules","blastn_stringtieRPF.snakefile")


# BLASTN stringtie RPF.fasta vs dd_Smed_v6.fasta
# ---------------------------------------------------------------------------------------
include: os.path.join(".","rules","final_novel_transcripts_prediction.snakefile")


# Predicting the fragemnted transcripts
# ---------------------------------------------------------------------------------------
include: os.path.join(".","rules","Fragmented.snakefile")
