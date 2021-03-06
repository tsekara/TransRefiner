smed_first_100=read.delim("~/owncloud.gwdg.de/Planarian_Simulation/Smed_cds_gmap_cov_0.8_idn_0.8/Extensions_first_100nts_deleted_intron_length_250.bed",header=F)
head(smed_first_100)
colnames(smed_first_100)=c("Sacffold","Extn_start","Extn_End","Counts", "Genome_strand" ,"Gene","Chopped_start","Chopped_end","Transcriptome_strand")
head(smed_first_100)
library(dplyr)
first_50_grouped=group_by(smed_first_100,Gene)
head(first_50_grouped)
head(first_50_grouped)
head(as.data.frame(summarize(sa,start=min(Extn_start),max=max(Extn_End))))
head(as.data.frame(summarize(first_50_grouped,start=min(Extn_start),max=max(Extn_End))))
first_100=as.data.frame(summarize(first_50_grouped,start=min(Extn_start),max=max(Extn_End)))
smed_orig_start_end=read.delim("/home/tsekara/owncloud.gwdg.de/Planarian_Simulation/Smed_cds_gmap_cov_0.8_idn_0.8/smed_genes_original_start_end.txt",header=F)
head(smed_orig_start_end)
colnames(smed_orig_start_end)=c("Gene","Orig_start","Orig_end")
head(merge(first_100,smed_orig_start_end,by="Gene"))
first_100=merge(first_100,smed_orig_start_end,by="Gene")
colnames(first_100)[2]
colnames(first_100)[2]="Extn_start"
colnames(first_100)[3]="Extn_end"
head(first_100)
head(first_50_grouped)
head(as.data.frame(summarize(sa,Extn_start=min(Extn_start),Extn_end=max(Extn_End),Genome_strand=Genome_strand)))
ls(=)
ls()
heaD(smed_first_100)
head(smed_first_100)
head(as.data.frame(summarize(smed_first_100,Extn_start=min(Extn_start),Extn_end=max(Extn_End),Genome_strand=Genome_strand)))
head(as.data.frame(summarize(smed_first_100,Extn_start=min(Extn_start),Extn_end=max(Extn_End),Genome_strand)))
rm(list=ls()=
rm(list=ls())
smed_orig_start_end=read.delim("/home/tsekara/owncloud.gwdg.de/Planarian_Simulation/Smed_cds_gmap_cov_0.8_idn_0.8/smed_genes_original_start_end.txt",header=F)
head(smed_orig_start_end)
smed_first_100=read.delim("~/owncloud.gwdg.de/Planarian_Simulation/Smed_cds_gmap_cov_0.8_idn_0.8/Extensions_first_100nts_deleted_intron_length_250.bed",header=F)
colnames(smed_first_100)=c("Sacffold","Extn_start","Extn_End","Counts", "Genome_strand" ,"Gene","Chopped_start","Chopped_end","Transcriptome_strand")
library(dplyr)
colnames(smed_first_100)=c("Gene","Orig_start","Orig_end","Strand")
head(smed_first_100)
colnames(smed_first_100)=c("Sacffold","Extn_start","Extn_End","Counts", "Genome_strand" ,"Gene","Chopped_start","Chopped_end","Transcriptome_strand")
colnames(smed_first_100)=c("Gene","Orig_start","Orig_end","Strand")
first_100_grouped=group_by(smed_first_100,Gene)
head(first_100_grouped)
head(smed_first_100)
colnames(smed_first_100)=c("Sacffold","Extn_start","Extn_End","Counts", "Genome_strand" ,"Gene","Chopped_start","Chopped_end","Transcriptome_strand")
head(smed_first_100)
first_100_grouped=group_by(smed_first_100,Gene)
head(first_100_grouped)
head(first_100_grouped)
dim(as.data.frame(summarize(first_100_grouped,Extn_start=min(Extn_start),Extn_end=max(Extn_End))))
head(as.data.frame(summarize(first_100_grouped,Extn_start=min(Extn_start),Extn_end=max(Extn_End)))))
head(as.data.frame(summarize(first_100_grouped,Extn_start=min(Extn_start),Extn_end=max(Extn_End))))
head(as.data.frame(summarize(first_100_grouped,Extn_start=min(Extn_start),Extn_end=max(Extn_End),Strand=as.character(Genome_strand))))
head(as.data.frame(summarize(first_100_grouped,Extn_start=min(Extn_start),Extn_end=max(Extn_End),Strand=as.character(Genome_strand))))
head(as.data.frame(summarize(first_100_grouped,Extn_start=min(Extn_start),Extn_end=max(Extn_End),Strand=as.integer(Genome_strand))))
head(as.data.frame(summarize(first_100_grouped,Extn_start=min(Extn_start),Extn_end=max(Extn_End))))
first_100=as.data.frame(summarize(first_100_grouped,Extn_start=min(Extn_start),Extn_end=max(Extn_End)))
head(first_100)
head(merge(first_100,smed_orig_start_end,by="Gene"))
head(smed_orig_start_end)
colnames(smed_orig_start_end)=c("Gene","Extn_start","Extn_end","Genome_strand")
head(smed_orig_start_end)
head(merge(first_100,smed_orig_start_end,by="Gene"))
colnames(smed_orig_start_end)=c("Gene","Orig_start","Orig_end","Genome_strand")
head(merge(first_100,smed_orig_start_end,by="Gene"))
first_100=merge(first_100,smed_orig_start_end,by="Gene"))
first_100=merge(first_100,smed_orig_start_end,by="Gene")
dim(first_100[(first_100$Genome_strand == "+" & first_100$Extn_start <  first_100$Orig_start | first_100$Genome_strand == "-" & first_100$Extn_end > first_100$Orig_end ), ])
write.table(first_100,"Extensions_and_orignal_first_100nts_Smed_genes.txt",sep="\t",quote=F,rownames=F)
write.table(first_100,"Extensions_and_orignal_first_100nts_Smed_genes.txt",sep="\t",quote=F,row.names=F)
smed_last_100=read.delim("~/owncloud.gwdg.de/Planarian_Simulation/Smed_cds_gmap_cov_0.8_idn_0.8/Last_100nts_deleted/Extensions_last_100nts_deleted_intron_length_200.bed",header=F)
colnames(smed_last_100)=c("Sacffold","Extn_start","Extn_End","Counts", "Genome_strand" ,"Gene","Chopped_start","Chopped_end","Transcriptome_strand")
head(smed_last_100)
last_100_grouped=group_by(smed_last_100)
head(as.data.frame(summarize(last_100_grouped,start=min(Extn_start),Extn_End=max(Extn_End))))
head(last_100_grouped)
head(as.data.frame(summarize(last_100_grouped,start=min(Extn_start),Extn_End=max(Extn_End))))
last_100_grouped=group_by(smed_last_100,Gene)
head(as.data.frame(summarize(last_100_grouped,start=min(Extn_start),Extn_End=max(Extn_End))))
last_100_summarized=as.data.frame(summarize(last_100_grouped,start=min(Extn_start),Extn_End=max(Extn_End)))
last_100_summarized=merge(last_100_summarized,smed_orig_start_end,by="Gene")
head(last_100_summarized)
colnames(last_100_summarized)[2]="Extn_start"
head(last_100_summarized)
dim(last_100_summarized[(last_100_summarized$Strand == "-" & last_100_summarized$Extn_start <  last_100_summarized$Orig_start | last_100_summarized$Genome_strand == "+" & last_100_summarized$Extn_End > last_100_summarized$Orig_end ), ])
head(last_100_summarized)
dim(last_100_summarized[(last_100_summarized$Strand == "-" & last_100_summarized$Extn_start <  last_100_summarized$Orig_start ), ])
dim(last_100_summarized[(last_100_summarized$Strand == "-" & last_100_summarized$Extn_start <  last_100_summarized$Orig_start ), ])
dim(last_100_summarized)
dim(last_100_summarized[(last_100_summarized$Strand == "-" & ), ])
dim(last_100_summarized[(last_100_summarized$Strand == "-"  ), ])
dim(last_100_summarized[(last_100_summarized$Genome_strand == "-" & last_100_summarized$Extn_start <  last_100_summarized$Orig_start | last_100_summarized$Genome_strand == "+" & last_100_summarized$Extn_End > last_100_summarized$Orig_end ), ])
dim(last_100_summarized)
last_100_summarized=last_100_summarized[(last_100_summarized$Genome_strand == "-" & last_100_summarized$Extn_start <  last_100_summarized$Orig_start | last_100_summarized$Genome_strand == "+" & last_100_summarized$Extn_End > last_100_summarized$Orig_end ), ]
head(last_100_summarized)
write.table(last_100_summarized,"/home/tsekara/owncloud.gwdg.de/Planarian_Simulation/Smed_cds_gmap_cov_0.8_idn_0.8/Last_100nts_deleted/Extensions_vs_original_last_100_smed_genes",sep="\t",row.names = F, quote=F)
ls()
head(first_100)
head(smed_first_100)
dim(smed_first_100)
dim(first_100)
head(merge(first_100,smed_first_100,by="Gene"))
dim(merge(first_100,smed_first_100,by="Gene"))
head(smed_orig_start_end)
head(smed_first_100)
smed_first_100_scaffold_gene$scaffold=smed_first_100$Sacffold
smed_scaffold_gene=()
smed_scaffold_gene=0
smed_scaffold_gene$scaffold=smed_first_100$Sacffold
smed_orig_start_end=read.delim("~/owncloud.gwdg.de/Planarian_Simulation/Smed_cds_gmap_cov_0.8_idn_0.8/smed_genes_original_start_end.txt",header =F)
head(smed_orig_start_end)
colnames(smed_orig_start_end)=c("Gene","Orig_start","Orig_end","Genome_strand","Scaffold")
head(smed_orig_start_end)
head(merge(first_100, smed_orig_start_end))
dim(merge(first_100, smed_orig_start_end))
first_100=merge(first_100, smed_orig_start_end)
write.table(first_100,"/home/tsekara/owncloud.gwdg.de/Planarian_Simulation/Smed_cds_gmap_cov_0.8_idn_0.8/First_100_deleted/Extensions_vs_original_first_100_smed_genes",sep="\t",row.names = F, quote=F)
head(first_100)
dim(merge(last_100_summarized, smed_orig_start_end))
head(last_100_summarized)
dim(last_100_summarized)
head(merge(last_100_summarized, smed_orig_start_end))
last_100_summarized=merge(last_100_summarized, smed_orig_start_end)
head(last_100_summarized)
write.table(last_100_summarized,"/home/tsekara/owncloud.gwdg.de/Planarian_Simulation/Smed_cds_gmap_cov_0.8_idn_0.8/First_100_deleted/Extensions_vs_original_first_100_smed_genes.bed",sep="\t",row.names = F, quote=F)
write.table(first_100,"/home/tsekara/owncloud.gwdg.de/Planarian_Simulation/Smed_cds_gmap_cov_0.8_idn_0.8/First_100_deleted/Extensions_vs_original_first_100_smed_genes",sep="\t",row.names = F, quote=F)
write.table(last_100_summarized,"/home/tsekara/owncloud.gwdg.de/Planarian_Simulation/Smed_cds_gmap_cov_0.8_idn_0.8/Last_100nts_deleted/Extensions_vs_Original_last_100_deleted_smed_genes.bed",sep="\t",row.names = F, quote=F)
savehistory("~/owncloud.gwdg.de/Planarian_Simulation/Smed_cds_gmap_cov_0.8_idn_0.8/Planarians_simualtion.R")
