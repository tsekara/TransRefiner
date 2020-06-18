# First, we load the libraries we need to carry out the analysis and generate the plots.
library("DESeq2")
library("scales")
library("RColorBrewer")
library("gplots")
library("calibrate")
# We assign the name of the conditions to variables. This way we can use the same script for other analyses by doing some minor changes.
s1<-"P6Cont"
s2<-"MafbcKO"


# Reading the metadata file:
meta <- read.delim(file = paste("meta_", s1, "_", s2, ".txt", sep=""), header = F)

# Formatting the table we have just read. This will tell DeSeq what each column means.
sampleTable<-data.frame(sampleName=meta$V1, fileName=meta$V2, condition=meta$V3)
options(stringsAsFactors = FALSE)

# Reading the data. Note that so far we have just loaded the metadata file. 
# The next command uses the data stored in our metadata file to load the actual counts for all samples.
# Importantly, we are establishing the design of the experiment; that is, we want to compare by condition.
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = ".", design= ~ condition)

# We supply the levels in order to make sure that the controls are the first elements.
colData(ddsHTSeq)$condition <- factor(colData(ddsHTSeq)$condition, levels=c(s1,s2))
colData(ddsHTSeq)$condition <- relevel(colData(ddsHTSeq)$condition, s1)

# To runDESeq, we just need one line:
ddsHTSeq <- DESeq(ddsHTSeq)

#########################
###   Quality Check   ###
#########################

### To begin with, we will mostly be using this DESeq object to do some quality checks

# DESeq works with the count data per gene. The data is not normalized between different libraries.
# We perform a Variance stabilizing transformation. This will facilitate comparisons between samples:
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)


### Heatmap of transformed values across samples.

# We choose the colours for the heatmap.
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

# To generate the heatmap we will use the genes with the greatest count across all samples (mean).
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:100]

# We open a pdf file
pdf(paste("plots/HM-vsd_",s2, ".pdf", sep=""))

# To do the heatmap:
heatmap.2(assay(vsd)[select,], col = hmcol, Rowv = FALSE, trace="none", margin=c(10, 6))

# Close the pdf
dev.off()


### Plotting the Euclidean distances between the samples, with transformed data.

# Calculate the distance across samples. This operation will be done on a per-column basis, so we transpose the matrix.
distsRL <- dist(t(assay(vsd)))

# Formatting the results of the previous step
mat <- as.matrix(distsRL)

# Naming the rows (so they are properly labeled in the plot).
rownames(mat) <- colnames(mat)

# Open pdf file
pdf(paste("plots/samp-distance_vsd_",s2, ".pdf", sep=""))

# Generate heatmap.
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))

# Close pdf file
dev.off()


### Principal component plot

pdf(paste("plots/all_pca_",s2, ".pdf", sep=""))
print(plotPCA(vsd, intgroup="conditions", ntop=nrow(assay(vsd))))
dev.off()


### Correlation across sample  log+1

# We will do the correlation of counts across samples. 

# The DESeq object has many methods associated with it. This are different ways the data can be accessed.
# The next line allows to retrieve a matrix with the counts of all genes. This matriz has one row per gene and one column per sample.
counts<-counts(dds)

# We will use the log of the counts because values have a very wide range. There will probably be many 0's, to avoid infinite values we will add 1 to all counts.
# We can declare functions in the following way:
f<-function(x){ log(x+1)}

# To apply the function to all elements of the matrix
trans_log<-matrix(sapply(counts, f), nrow=nrow(counts), ncol=ncol(counts), byrow=F)

# Open the pdf file
pdf(paste("plots/log-1_correlation_across_samples_",s2, ".pdf", sep=""))

# Change the parameters to enable more than one plot to be printed in the same sheet.
# With mfrow we set the number of columns and rows to equal the number of samples. The rest are graphical parameters.  
par(mfrow=c(ncol(trans_log),ncol(trans_log)), oma=c(2,2,6,2), mar=rep(0.5, 4))

# Use a loop to compare all samples in a pairwise manner.
# The outer cycle will be executed N-sample of times. Variable "i" will take values from 1 to N-sample
for(i in 1:ncol(trans_log)){
  # The inner cycle behaves as the outer cycle. The result is that we will carry out the following instructions a total of N*N times.
  for(j in 1:ncol(trans_log)){
    # We do not want to compare a sample to itselt. We use this conditional clause to establish the behavior when variables i and j are equal.
    # In the diagonal we will print the name of the sample.	
    if(i==j){
      # We "draw" an empty plot
      plot(2, col="white",xaxt='n',yaxt='n',pch='',ylab='',xlab='')
      # We write on it the name of the sample.
      text(1,2,colnames(counts)[i], cex=2, col="purple")
    }
    # In the "lower" side of the diagonal we will print text as well, in this case, the value of the Spearman's rho.
    # That is, when the row number is smaller than the column number
    if(j<i){
      # We calculate the correlation of sample i and j with the spearman method.
      spearm<-cor.test(trans_log[,i], trans_log[,j], method="spearman")
      # We draw an empty plot
      plot(2,2, col="white", xaxt='n',yaxt='n',pch='',ylab='',xlab='')
      # I print the value of the rho, with two significance digits.
      text(2,2, signif(spearm$estimate[[1]], 2), cex=1.5)
    }
    # In the "upper" side of the diagonal we draw a scatterplot of the values of one sample in the y axis,
    # in the x axis we draw the corresponding values of the other sample being considered.
    if(j>i){
      # We plot the values for all genes
      plot(trans_log[,i], trans_log[,j], pch=16,  cex=0.7, col=alpha("black", 0.5), yaxt="n", xaxt="n")
      # We draw the regression line
      abline(lm(trans_log[,j]~trans_log[,i]), col="seagreen1")
    }
    if(i==1 & i!=j){
      # We draw the x axis on the first row of plots.
      axis(3)
    }
    if(j==ncol(trans_log) & i!=j){
      # We draw the y axis for the first column of plots
      axis(4)
    }
  }
}
# We add a title for our "grid" of plots
mtext("Correlation across samples",  NORTH<-3,line=3, adj=0.5, cex=1.2, outer=TRUE)
# Close the pdf
dev.off()



#############################################################################################


### Differential Gene expression analysis

# We have already run DeSEq. We can extract the results with the following line:
res <- results(ddsHTSeq)

# We filter out rowsthat have an "NA" in the adjusted p-value column
de.filtered = res[!is.na(res$padj), ]

# Running biomart to fetch gene names
library("biomaRt")
ensembl = useMart ("ensembl")

ensembl = useMart ("ensembl", dataset="mmusculus_gene_ensembl")
results <- getBM (attributes = c("ensembl_gene_id", "external_gene_id"),
                  filters = "ensembl_gene_id",
                  values = rownames (de.filtered),
                  mart = ensembl)

idx <- match (rownames (de.filtered), results$ensembl_gene_id)
de.filtered$gene_name <- results[idx, "external_gene_id"]


### Tables
# We can already print the table with the results
write.table (res, file = paste (s1, "-", s2, "_all_genes.txt", sep = ""), quote = F, row.names = F, col.names = T, sep = "\t", na = "")

# We can subset out results to include only significantly, differentially expressed genes.
significant.genes=de.filtered[(de.filtered$padj < .05),]

# Se can also use a threshold for a minimum fold change:
up<- significant.genes[(significant.genes$log2FoldChange > 0),]
down<- significant.genes[(significant.genes$log2FoldChange < 0),]

# We write the filtered results to a plain text file. These will be "tab" delimited tables:
write.table(up,file=paste("up_", s1,"-",s2,".genes.txt", sep=""), quote=F, row.names=T, col.names=T, sep = "\t")
write.table(down,file=paste("down_", s1,"-",s2,".genes.txt", sep=""), quote=F, row.names=T, col.names=T, sep = "\t")


### MA Plot 

color_dots <- "gray32"
color_de <- "pink"
color_marked <- "blue"

# Open pdf
pdf(paste("plots/", s1, "-", s2, "_ma_plot.pdf", sep=""))

# We use the MA plot function from DESeq2
plotMA(dseq_result, alpha=0.05)

# Close pdf
dev.off()


##### Volcano plot

# Open pdf
pdf(paste("plots/", s1, "-", s2, "_volcano.pdf", sep=""))

# We plot
plot(de.filtered$log2Fold,-log10(de.filtered$pval), pch=16, cex=0.6, col = ifelse((de.filtered$padj<=0.05 & abs(de.filtered$log2Fold)>=1), alpha("violetred1", 0.5), alpha("gray50", 0.5)), ylab="-log10 pvalue", xlab="log2 fold change", main=paste("Comparison of samples ", s2, "vs", s1, sep=" "))

# We can include lines that will help guide us.
abline(b=0, v=1, col=alpha ("black", 0.9))
abline(b=0, v=-1, col=alpha ("black", 0.9))

# We can add a legend to make the plot easier to understand.
legend("topright", "Genes with a log2Fold change > 1 and padj <0.05", x.intersp=0.3, pch=20, col=alpha(color_de, 0.5), cex=0.6, horiz=TRUE)

# Close pdf
dev.off()


# We don't always need to output the plot to a pdf, we can visialize it from R. 
plot(de.filtered$log2Fold,-log10(de.filtered$pval), pch=16, cex=0.6, col = ifelse((de.filtered$padj<=0.05 & abs(de.filtered$log2Fold)>=1), "violetred1", "gray50"), ylab="-log10 pvalue", xlab="log2 fold change", main=paste("Comparison of samples ", slit, "vs", control, sep=" "))
abline(b=0, v=1, col="darkorchid1", lty=2)
abline(b=0, v=-1, col="darkorchid1", lty=2)
legend("topright", "DE genes", x.intersp=0.3, pch=20, col=alpha("violetred1", 0.5), cex=0.6, horiz=TRUE)

# We can highlight some gene of interest:
textxy(de.filtered["ENSMUSG00000074622",2], -log10(de.filtered["ENSMUSG00000074622",5]), "Mafb", cex=0.9, offset=0.6, col="steelblue4")
points(de.filtered["ENSMUSG00000074622",2], -log10(de.filtered["ENSMUSG00000074622",5]), col="steelblue4", cex=0.6, pch=16)

# Saving the session
save.image("de.RData")