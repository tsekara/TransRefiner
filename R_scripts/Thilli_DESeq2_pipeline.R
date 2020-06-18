##################################################################### Loading Libraries #############################################################################################
library("DESeq2")
library("scales")
library("RColorBrewer")
library("gplots")
library("calibrate")
library("pheatmap")
library("pca3d")
library("geneplotter")

##################################################################### Loading scripts ##########################################################################
source("/home/tsekara/R_Scripts/Scatterplots_tiff.R")
source("/home/tsekara/R_Scripts/Scatterplots_pdf.R")

################################################################ CLEAR the HISTORY #####################################################################################

########################################################################################################################################################################################
########################################################## Converting sam2profile counts to count-table ###############################################################################
########################################################################################################################################################################################

count.df <- read.csv("");
count.mat <- xtabs(Reads ~ Target + Experiment, data = count.df);
count.mat <- count.mat[order(rowSums(count.mat), decreasing = TRUE),];
write.table(count.mat,"");

########################################################################################################################################################################################
########################################################## Reading the count data genrated in provious tep into R ######################################################################
########################################################################################################################################################################################

data = read.delim("xxxx.txt", row.names=1)

                        ############################################## Scatterplots(Tiff) #####################################

scatterplots_tiff(data)

                       ########################################### Scatterplots(PDF) #########################################
scatterplots_pdf(data)

                      ############################################ generating matrix ##########################################
        
rnaseqMatrix = as.matrix(round(data))


#rnaseqMatrix = rnaseqMatrix[rowSums(rnaseqMatrix)>=10,]

########################################################################## Defining Conditions #######################################################################################

conditions = data.frame(conditions=factor(c(rep("GFP", 3), rep("H2B", 3))))
rownames(conditions) = colnames(rnaseqMatrix)
ddsFullCountTable <- DESeqDataSetFromMatrix(countData = rnaseqMatrix,colData = conditions,design = ~ conditions)
dds = DESeq(ddsFullCountTable)

############################################################ Dispersion plot and Size factor estimation #####################################################

sizeFactors(dds)
plotDispEsts(dds)

############################################################### How Normlaization works ##################################################################################

GeneCounts <- counts(dds)
idx.nz <- apply(GeneCounts, 1, function(x) { all(x > 0)})  

multidensity( counts(dds, normalized = T)[idx.nz ,],xlab="mean counts", xlim=c(0, 1000))

########################################################################################################################################################################################
############################################################################# QUALITY CONTROL #######################################################################################
########################################################################################################################################################################################

################################################################################## VSD transformation ########################################################################        

vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vsd_counts=assay(vsd)

      
########################################################################################################################################################################################
########################################################## Sample-to-Sample distance ######################################################################
########################################################################################################################################################################################

                                   ###################### Caluculating the distance ################

Sample_dist=dist(t(log(vsd_counts)))
Sample_dist_mat=as.matrix(Sample_dist)                 

                                   ###################### Computing Hierachicahl Clustering #########

hclust_Sample=hclust(Sample_dist)    

                                    ###################### Sample-Sample distance heatmap ###########

pheatmap(Sample_dist_mat,Colv=as.dendrogram(hclust_Sample),Rowv = as.dendrogram(hclust_Sample),trace="none",margin=c(10,10),col=heat.colors(256),annotation = conditions)

                                    ####################### Correlation Plot #########################################

pheatmap(cor(vsd_counts),Colv=as.dendrogram(hclust_Sample),Rowv = as.dendrogram(hclust_Sample),trace="none",margin=c(10,10),annotation = conditions,border_color = NA)

########################################################################################################################################################################################
################################################################################################### PCA 2D & 3D ######################################################################
########################################################################################################################################################################################


####################################################### 2D PCA Plot ##########################################################################
print(plotPCA(vsd, intgroup="conditions", ntop=nrow(assay(vsd))))

                                             ####################### or ##############################
pr=prcomp(t(vsd_counts))
plot(pr$x,col="red",main="PC-plot",xlim=c(-22,15))
text(pr$x[,1],pr$x[,2],labels=colnames(vsd_counts),cex=.7)

###################################################### 3D PCA plot ##########################################################################
pr=prcomp(t(vsd_counts))
pca3d(pr,1:3,show.labels = T,show.axes = T)

########################################################################################################################################################################################
################################################################################## Differential Gene Expression ######################################################################
########################################################################################################################################################################################
# res = results(dds,independentFiltering = FALSE)
res = results(dds,independentFiltering = TRUE)

baseMeanA <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$condition == "GFP"])
baseMeanB <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$condition == "H2B"])
res = cbind(GFP=baseMeanA,H2B=baseMeanB, as.data.frame(res))
res = cbind(id=rownames(res), as.data.frame(res))

write.table(as.data.frame(res[order(res$padj),]), file='DEG_without_annotations.txt', sep='\t', quote=FALSE, row.names=F)



sig  <- res[which(res$padj  < 0.05),]
sig.deseq  <- rownames(sig)

################################################################## Calling annotation ###############################################################################################
tr5_frameDP=read.delim("/home/tsekara/Transcriptome_Annotations/Trinity5_FrameDP_corrected/Annotations/tr5_frameDPcorrected_annotations_v1.txt")
dresden_annotations=read.delim("/home/tsekara/Transcriptome_Annotations/Dresden_transcriptome/dd_annot_v4.txt")

################################################################# Merging trinity5_frameDP Annotations and expression data ###########################################################################

tr5_DEG_annotations=merge(res,tr5_frameDP,by.x="id",by.y="Transcript_ID")
write.table(as.data.frame(tr5_DEG_annotations[order(tr5_DEG_annotations$padj),]), file='tr5_DEG_annotations.txt', sep='\t', quote=FALSE, row.names=F)

################################################################# Merging dresden Annotations and expression data ###########################################################################
dd_DEG_annotations=merge(res,dresden_annotations,by.x="id",by.y="Transcript_ID")
write.table(as.data.frame(dd_DEG_annotations[order(dd_DEG_annotations$padj),]), file='dd_DEG_annotations.txt', sep='\t', quote=FALSE, row.names=F)

########################################################################################################################################################################################
################################################################################## Volcano plot and MA plot ######################################################################
########################################################################################################################################################################################

################################################ Volcano Plot #####################################################
plot(res$log2FoldChange,-log10(res$padj),pch=16,ylab="-Log10_padj", xlab="Log2_fold_change",main="Volcano Plot",family="sans", font.axis=2, cex.axis=2, cex.lab=1.25, font.lab=2,col=ifelse(res$padj<=0.05 & res$log2FoldChange >1, "red",ifelse(res$padj<=0.05 & -res$log2FoldChange >1 ,"blue","black")))
legend("topright", c("LFC > 1 & padj <0.05","LFC < -1 & padj <0.05"), x.intersp=1, pch=20, col=c("red","blue"), cex=1.5, horiz=F)

                                                              OR

tiff("Volacno.tiff", height = 12, width = 17, units = 'cm',compression = "lzw", res = 400)
plot(res$log2FoldChange,-log10(res$padj),pch=16,ylab="-Log10_padj", xlab="Log2_fold_change",main="Volcano Plot",family="sans", font.axis=2, cex.axis=2, cex.lab=1.25, font.lab=2,col=ifelse(res$padj<=0.05 & res$log2FoldChange >1, "red",ifelse(res$padj<=0.05 & -res$log2FoldChange >1 ,"blue","black")))
legend("topright", c("LFC > 1 & padj <0.05","LFC < -1 & padj <0.05"), x.intersp=1, pch=20, col=c("red","blue"), cex=1.5, horiz=F)
dev.off()
dev.off()

plot(res$log2FoldChange,-log10(res$padj),pch=16,ylab="-log10_padj", xlab="log2_fold_change")
points(sig$log2FoldChange,-log10(sig$padj),col="red",pch=16)
abline(b=0, v=.5, col=alpha ("black", 0.9))
abline(b=0, v=-.5, col=alpha ("red", 0.9))

################################################ MA Plot #####################################################
res$padj[is.na(res$padj)]  <- 1


plot(log2(res$baseMean), res$log2FoldChange,pch=16)
plot(log2(res$baseMean), res$log2FoldChange,pch=16, col=ifelse(res$padj<=0.05 & res$log2FoldChange >1, "red",ifelse(res$padj<=0.05 & -res$log2FoldChange >1 ,"blue","black")), xlab="Log2(res$baseMean)", ylab="Log2FoldChange", main="MA-plot",family="sans", font.axis=2, cex.axis=2, cex.lab=1.25, font.lab=2);
legend("topright", c("LFC > 1 & padj <0.05","LFC < -1 & padj <0.05"), x.intersp=1, pch=20, col=c("red","blue"), cex=1.5, horiz=F)


tiff("MA.tiff", height = 12, width = 17, units = 'cm',compression = "lzw", res = 400)
plot(log2(res$baseMean), res$log2FoldChange,pch=16, col=ifelse(res$padj<=0.05 & res$log2FoldChange >1, "red",ifelse(res$padj<=0.05 & -res$log2FoldChange >1 ,"blue","black")), xlab="Log2(res$baseMean)", ylab="Log2FoldChange", main="MA-plot",family="sans", font.axis=2, cex.axis=2, cex.lab=1.25, font.lab=2);
legend("topright", c("LFC > 1 & padj <0.05","LFC < -1 & padj <0.05"), x.intersp=1, pch=20, col=c("red","blue"), cex=1.5, horiz=F)
dev.off()
dev.off()

                                                        #    OR

                                      ############## MAplot with basemean_mRNA in X-axis and log2_foldchange_RPF in yaxis  ################
plot(log2(res$baseMean_mRNA), res$log2FoldChange,pch=16)
plot(log2(res$baseMean_mRNA), res$log2FoldChange,pch=16, col=ifelse(res$padj<=0.05 & res$log2FoldChange >1, "red",ifelse(res$padj<=0.05 & -res$log2FoldChange >1 ,"blue","black")), xlab="Log2(res$baseMean_mRNA)", ylab="Log2FoldChange", main="MA-plot",family="sans", font.axis=2, cex.axis=1.25, cex.lab=2, font.lab=2);

                                                                                    OR

plot(log2(res$baseMean), res$log2FoldChange, col=ifelse(res$padj<=0.05 ,"red","black"), xlab="log-Counts", ylab="Log2FoldChange", main="MA-plot");
abline(h=.5, col=alpha ("blue", 0.9))
abline(h=-.5, col=alpha ("green", 0.9))
