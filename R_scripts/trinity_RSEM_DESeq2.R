library(DESeq2)
######## Importing the raw read counts matrix estimated by RSEM #######################################

data = read.table("/home/tsekara/Condition_slit_dd_smed/./matrices_ohne_Condition_1/Tr_ohne_Condition_1.counts.matrix", header=T, row.names=1, com='')

######### Definig what columns has to make pairwise comparisions from raw matrix ###########

col_ordering = c(1,2,3,4,5,6)
rnaseqMatrix = data[,col_ordering]
rnaseqMatrix = round(data)
rnaseqMatrix = rnaseqMatrix[rowSums(rnaseqMatrix)>=2,]

conditions = data.frame(conditions=factor(c(rep("Control", 3), rep("DNAH11", 3))))

rownames(conditions) = colnames(rnaseqMatrix)

ddsFullCountTable <- DESeqDataSetFromMatrix(countData = rnaseqMatrix,colData = conditions,design = ~ conditions)
dds = DESeq(ddsFullCountTable)

res = results(dds,independentFiltering = FALSE)

baseMeanA <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$condition == "Control"])
baseMeanB <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$condition == "DNAH11"])

res = cbind(baseMeanA, baseMeanB, as.data.frame(res))

res = cbind(id=rownames(res), as.data.frame(res))

res$padj[is.na(res$padj)]  <- 1

write.table(as.data.frame(res[order(res$pvalue),]), file='DEG_with all_replicates.txt', sep='\t', quote=FALSE, row.names=F)

###################### Calling rnaseq_plot functions R script ###################################
source("/home/tsekara/softwares/trinityrnaseq/Analysis/DifferentialExpression/R/rnaseq_plot_funcs.R")

pdf("Counts.matrix.Control_vs_Condition.DESeq2.DE_results.MA_n_Volcano.pdf")

plot_MA_and_Volcano(log2(res$baseMean+1), res$log2FoldChange, res$padj)

dev.off()