dmel_cds_gene_start_end=read.delim("dmel_CDS_start_end_strand.txt",header=F)
colnames(dmel_cds_gene_start_end)=c("Gene","Orig_Start","Orig_end","Strand")
first_100=read.delim("../First_100nts_deleted/Cleaned_Extensions_first_100nts_deleted_intron_length_70.bed",header=F)
colnames(first_100)=c("Chromosome","Extn_Start","Extn_End","Counts","Genome_strand", "Gene","Transcript_Start","Transcript_End","Strand")
first_100_group=group_by(first_100,Gene)
first_100nts=as.data.frame(summarize(first_100_group,Extn_start=min(Extn_Start),Extn_end=max(Extn_End)))
dim(first_100nts)
first_100nts=merge(first_100nts,dmel_cds_gene_start_end,by="Gene")
Extn_greater_orig_first_100nts=first_100nts[(first_100nts$Strand == "+" & first_100nts$Extn_start <  first_100nts$Orig_Start
                                             |first_100nts$Strand == "-" & first_100nts$Extn_end > first_100nts$Orig_end ), ]
dim(Extn_greater_orig_first_100nts)
Extn_lesser.equal_orig_first_100nts=first_100nts[(first_100nts$Strand == "+" & first_100nts$Extn_start >=  first_100nts$Orig_Start
                                             |first_100nts$Strand == "-" & first_100nts$Extn_end <= first_100nts$Orig_end ), ]
dim(Extn_lesser.equal_orig_first_100nts)


