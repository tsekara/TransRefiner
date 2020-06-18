                      ############################################## DESeq2 ########################################################

                      ###################### importng the count reads from samttols idxcountsfrom all samples ######################
                      library(DESeq2)
                      counts.df=NULL               
                      experiment.names <- c("control_1","control_2","control_3","slit_1","slit_2","slit_3","wnt5_1","wnt5_2","wnt5_3");
                      for(expName in experiment.names)
                         {
                             reads.df <- read.delim(sprintf("%s_idxcounts.txt",expName),col.names = c("refName","refLength","mapped","unmapped"));
                             counts.df <- cbind(counts.df, reads.df$mapped);
                         }
                      rownames(counts.df) <- reads.df$refName;
                      colnames(counts.df) <- experiment.names;                      
                      ################################## loading the samplesheet ####################################
                      samplesheet=data.frame(row.names=colnames(counts.df),condition=c("control","control","control","slit","slit","slit","wnt5","wnt5","wnt5"))
                      ################################## Creating a Condition factor ########################################
                      condition=factor(c("control","control","control","slit","slit","slit","wnt5","wnt5","wnt5")) 
                      ################################# loading the DESEq2dataset #####################################
                      dds<-DESeqDataSetFromMatrix(countData= counts.df,colData=samplesheet,design=~condition)
                      ################################# Finding the differential expressed genes between control vs lhx1 ################################
                      dseq2_result=DESeq(dds)
                      ################################ creating the results ##########################################
                      result_control_slit=results(dseq2_result,contrast=c("condition","control","slit"))
                      result_control_wnt5=results(dseq2_result,contrast=c("condition","control","wnt5"))
                      ################################# Ordering the results based on Padj.values #############################################
                      results_padj=res[order(res$padj <.001),]
                      ################################# Ordering the results based on log2foldchange of 1 ####################################
                      res_lfc_up=res_pdj[order(res_pdj$log2FoldChange > 1),]
                      res_lfc_down=res_pdj[order(res_pdj$log2FoldChange < -1),]
                      ################################# rlog Transformation for H.Clustering, PCA, etc., ################################
                      rld=rlog(dds)
                      rlogMat=assay(rld)
                                                                    #### or ### ##### 
                                                                ###### which one is better?#######
                      vsd=varianceStabilizingTransformation(dds)
                      vstMat=assay(vsd)
                      ################################ Top 30 expressed geens ########################################
                      select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:30]
                      hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
                      heatmap.2(counts(dds,normalized=TRUE)[select,], col = hmcol,Rowv = FALSE, Colv = FALSE, scale="none",dendrogram="none", trace="none", margin=c(10,6))
                      ################################ PCA ##############################################################
                      plotPCA(rld,intgroup=c("condition"))
                      ############################## Distance matrix and Correlation plot ###################################
                      distsRL <- dist(t(assay(rld)))
                      mat <- as.matrix(distsRL)
                      rownames(mat) <- colnames(mat) <- with(colData(dds),paste(condition, type, sep=" : "))
                      rownames(mat) <- colnames(mat) <- with(colData(dds),paste(condition, sep=" : "))
                      hc <- hclust(distsRL)
                      heatmap.2(mat, Rowv=as.dendrogram(hc),symm=TRUE, trace="none", col = rev(hmcol), margin=c(13, 13))

