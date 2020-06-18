sp_GO=read.delim("C:/Users/tsekaran/Desktop/Owncloud/Jose_GO_Enrichment/Control_vs_random_sampling/max_blastx_hits_10/eval_10/Swisprot_oct_2016_GO_terms.txt",header=T)
head(sp_GO)
sp_GO$Entry.name=NULL
setwd("/Users/tsekaran/Desktop/Owncloud/Jose_GO_Enrichment/Control_vs_random_sampling/max_blastx_hits_10/eval_.001/")
blatsx=read.delim("dd_semd_v4_vs_uniprot_oct16_max_blastx_hits_10_eval_.001.txt",header=F)
colnames(blastx)
colnames(blatsx)
smed_dd_v4_vs_sp_id=read.delim("dd_semd_v4_vs_uniprot_IDs_oct16_max_blastx_hits_10_eval_.001")
head(smed_dd_v4_vs_sp_id)
head(merge(smed_dd_v4_vs_sp_id,sp_GO,by="Entry"))
smed_v4_vs_sp_GO_max_target_5=merge(smed_dd_v4_vs_sp_id,sp_GO,by="Entry")
smed_v4_vs_sp_GO_max_target_5$Entry=NULL
library(data.table)
res=setDT(smed_v4_vs_sp_GO_max_target_5)[, .(Gene.ontology.IDs= paste(unique(unlist(strsplit(Gene.ontology.IDs, "; "))), collapse="; ")), by = Transcript_Id ]
setDT(smed_v4_vs_sp_GO_max_target_5)[, .(Gene.ontology.IDs= paste(unique(unlist(strsplit(Gene.ontology.IDs, "; "))), collapse="; ")), by = Transcript_Id ]
head(merge(smed_dd_v4_vs_sp_id,sp_GO,by="Entry"))
head(smed_dd_v4_vs_sp_id)
head(smed_v4_vs_sp_GO_max_target_5)
write.table(smed_v4_vs_sp_GO_max_target_5,"a.txt",sep="\t",row.names=F,quote=F)
setDT(smed_v4_vs_sp_GO_max_target_5)[, .(Gene.ontology.IDs= paste(unique(unlist(strsplit(Gene.ontology.IDs, "; "))), collapse="; ")), by = Transcript_Id ]
setDT(smed_v4_vs_sp_GO_max_target_5)[, .(Gene.ontology.IDs= paste(unique(unlist(strsplit(Gene.ontology.IDs, ";"))), collapse="; ")), by = Transcript_Id ]
setDT(smed_v4_vs_sp_GO_max_target_5)[, .(Gene.ontology.IDs= paste(unique(unlist(strsplit(Gene.ontology.IDs, "; "))), collapse="; ")), by = Transcript_Id ]
smed_v4_vs_sp_GO_max_target_5=as.character(smed_v4_vs_sp_GO_max_target_5$Gene.ontology.IDs)
res=setDT(smed_v4_vs_sp_GO_max_target_5)[, .(Gene.ontology.IDs= paste(unique(unlist(strsplit(Gene.ontology.IDs, "; "))), collapse="; ")), by = Transcript_Id ]
colnames(smed_v4_vs_sp_GO_max_target_5)
head(smed_v4_vs_sp_GO_max_target_5)
blastx_sp=read.delim("a.txt")
blastx_sp$Gene.ontology.IDs=as.character(blastx_sp$Gene.ontology.IDs)
res=setDT(blastx_sp)[, .(Gene.ontology.IDs= paste(unique(unlist(strsplit(Gene.ontology.IDs, "; "))), collapse="; ")), by = Transcript_Id ]
write.table(res,"Summarized_smed_v4_vs_unique_uniprot_SP_GO_Ids_max_target_10_eval_.001.txt",sep="\t",row.names = F,quote = F)
library(splitstackshape)
install.packages("splitstackshape")
df2 <- cSplit(res, "Gene.ontology.IDs", sep = ";", direction = "long")
library(splitstackshape)
df2 <- cSplit(res, "Gene.ontology.IDs", sep = ";", direction = "long")
write.table(df2, "Bingo_Smed_v4_swisprot_max_target_10_eval_.001.txt",sep=" = ",quote=F,row.names=F)
save.image("Bingo_custom_db_generation_eval_.001_max_target_10.RData")
savehistory("Bingo_custom_db_generation_eval_.001_max_target_10.R")
