
dmel_cds_gene_start_end=read.delim("dmel_CDS_start_end_strand.txt",header=F)
colnames(dmel_cds_gene_start_end)=c("Gene","Orig_Start","Orig_end","Strand")
last_100=read.delim("Cleaned_Extensions_last_100nts_deleted_intron_length_70.bed",header=F)
colnames(last_100)=c("Chromosome","Extn_Start","Extn_End","Counts","Genome_strand", "Gene","Transcript_Start","Transcript_End","Strand")
last_100_group=group_by(last_100,Gene)
last_100nts=as.data.frame(summarize(last_100_group,Extn_start=min(Extn_Start),Extn_end=max(Extn_End)))
length(unique(last_100$Gene))
Total_last100_extension=length(unique(last_100$Gene)) 
last_100nts=merge(last_100nts,dmel_cds_gene_start_end,by="Gene")
last_100nts_extn_greater_orig=last_100nts[(last_100nts$Strand == "-" & last_100nts$Extn_start <  last_100nts$Orig_Start | last_100nts$Strand == "+" & last_100nts$Extn_end > last_100nts$Orig_end ), ]
dim(last_100nts_extn_greater_orig)
extensions_greater_than_orig_last_100nts=dim(last_100nts_extn_greater_orig)[1]
last_100nts_extn_less.equals_orig=last_100nts[(last_100nts$Strand == "-" & last_100nts$Extn_start >=  last_100nts$Orig_Start | last_100nts$Strand == "+" & last_100nts$Extn_end <= last_100nts$Orig_end ), ]
dim(last_100nts_extn_less.equals_orig)
extensions_greater_than_orig_last_100nts=dim(last_100nts_extn_less.equals_orig)[1]