############## Creating a Gene_to_trans_Map #######################
gene_trans_map=read.delim("/path/to/transcript-ids", header=F)
library(stringr)
head(str_extract(gene_trans_map$V1, "^dd_Smed_v6_[0-9]+_[0-9]+"))
gene_trans_map$V2=str_extract(gene_trans$V1, "^dd_Smed_v6_[0-9]+_[0-9]+")
gene_trans_map_smed_dd_v6_complete=gene_trans[c(2,1)]
write.table(gene_trans_map_smed_dd_v6_complete,"gene_trans_map_smed_dd_v6_complete.txt",sep="\t",row.names=F,quote=F,col.names = F)