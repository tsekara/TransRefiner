###### Loading required libraries ########################
library("DESeq2")
##### Calling in the sam2prfile counts ##################
reads.3col.df <- read.csv("target_counts_same2profile_output.csv");
experiment.counts <- xtabs(Reads ~ Target + Experiment, data = reads.3col.df);
(colnames(experiment.counts) <- sub("^Sample","",colnames(experiment.counts)));
######## Creating Experiment matrix ######################
experiment.tech.conditions <- sub("^(.*?)[ab]?$","\\1",colnames(experiment.counts));
experiment.tech.counts <- NULL;
for(x in unique(experiment.tech.conditions)){
  if(length(grep(x,colnames(experiment.counts))) > 1){
    experiment.tech.counts <-
      cbind(experiment.tech.counts,
            rowSums(experiment.counts[,grep(x,colnames(experiment.counts))]));
  } else {
    experiment.tech.counts <- cbind(experiment.tech.counts, experiment.counts[,x]);
  }
  colnames(experiment.tech.counts)[dim(experiment.tech.counts)[2]] <- sub("Sample","",x);
}
####### Samplesheet creation ####################    
samplesheet=data.frame(row.names=colnames(counts.df),condition=c("control","control","control","slit","slit","slit","wnt5","wnt5"))
######Creating a Condition factor #########
condition=factor(c("control","control","control","slit","slit","slit","wnt5","wnt5")) 
#### Loading the DESEq2dataset ###################
dds<-DESeqDataSetFromMatrix(countData=experiment.tech.counts,colData=samplesheet,design=~condition)
dds$condition=relevel(dds$condition,“control”)
######### Finding the differential expressed transcripts ########################
dseq2_result=DESeq(dds)
######### Creating the results ############################ 
result_control_slit=results(dseq2_result,contrast=c("condition","control","slit"),independentFiltering=FALSE)
result_control_wnt5=results(dseq2_result,contrast=c("condition","control","wnt5"),independentFiltering=FALSE)
######### Ordering the results based on Padj.values ########################## 
results_padj=res[order(res$padj <.05),]
######### Ordering the results based on log2foldchange of 1 ###################### 
res_lfc_up=res_pdj[order(res_pdj$log2FoldChange > 1),]
res_lfc_down=res_pdj[order(res_pdj$log2FoldChange < -1),]
######## rlog Transformation for H.Clustering, PCA, etc., ########################
rld=rlog(dds)
rlogMat=assay(rld)
### or #############
vsd=varianceStabilizingTransformation(dds)
vstMat=assay(vsd)
########## Top 30 expressed genes ###############################
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:30]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(counts(dds,normalized=TRUE)[select,], col = hmcol,Rowv = FALSE, Colv = FALSE,scale="none",dendrogram="none", trace="none", margin=c(10,6))
########################## PCA #############################
plotPCA(rld,intgroup=c("condition"))
###### Distance matrix and Correlation plot ####################################
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds),paste(condition, type, sep=" : "))
rownames(mat) <- colnames(mat) <- with(colData(dds),paste(condition, sep=" : "))
hc <- hclust(distsRL)
heatmap.2(mat, Rowv=as.dendrogram(hc),symm=TRUE, trace="none", col = rev(hmcol), margin=c(13, 13))