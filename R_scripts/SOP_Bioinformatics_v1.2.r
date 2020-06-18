###
# script for standard Illumina microarray analysis
###

# author Andrea Hofmann
# March 2011
# version 1.0

# Start MA analysis

# Only use this when installing R the first time
source("http://bioconductor.org/biocLite.R");
biocLite();
install.packages( c("RColorBrewer", "amap", "pcurve", "survival","vegan","ade4","scatterplot3d"));
source("http://www.bioconductor.org/biocLite.R");
biocLite("made4");

# start analysis
setwd("F:/Interne Projekte/SOP")

#select right annotation file from: 
#HumanHT-12_V4, HumanHT-12_V3, HumanMI_V1, HumanMI_V2, HumanWG-6_V2, HumanWG-6_V3, MouseMI_V1, MouseMI_V2, MouseWG-6_V1, MouseWG-6_V2
annotation <- read.delim("AnnotationIllumina/MouseWG-6_V2.txt",row.names="Array_Address_Id",dec=",")

load("Scripts_update_v1_2.RData")
# Read in data
# check for comma or point as dec
raw.Data <- read.delim("TestData.txt", row.names = 1, dec = ".") 
head(raw.Data)

# raw expression data
raw.expression <- raw.Data[,seq(1,dim(raw.Data)[2],2)]
dim(raw.expression)
# raw detection p value
raw.calls <- raw.Data[,seq(2,dim(raw.Data)[2],2)]
dim(raw.calls)

# Sample sheet needs: Array_ID	Lab_ID	Replicate	Group 
samplesheet <- read.delim("SampleSheet.txt")
# compare samplesheet to arraynames  
# => should be integer(0) = 0 otherwise return position of wrong annotation
compare.samplesheet(colnames(raw.expression),samplesheet$"Array_ID")

colnames(raw.expression) <- samplesheet$"Group"
colnames(raw.calls) <- samplesheet$"Group"

# QUALITY CONTROL
# boxplot
boxplot(raw.expression,log="y", pch='.',col=rainbow(ncol(raw.expression))) 
# scatterplots raw: create folder QC_raw  
create.scatterplot(raw.expression,"QC_raw")
# present call rate
present <- vector()
for(i in 1:ncol(raw.calls)){
  present[i] <- sum(raw.calls[,i]<0.05)/nrow(raw.calls)*100
}
writetable(cbind(samplesheet,present),"present",p_row=F)

# normalize data with quantiles normalization                                      
library (limma)
data.quantiles<- normalizeBetweenArrays (as.matrix (raw.expression), method = "quantile")
colnames (data.quantiles) <- colnames (raw.expression)
rownames (data.quantiles) <- rownames (raw.expression)

# create folder QC_norm
create.scatterplot(data.quantiles,"QC_norm")

# correlation matrix
correlation.heatmap(data.quantiles)

# variable genes
vargenes_0.5_10 <- getVariableGenes(data.quantiles,0.5,10)
sum(vargenes_0.5_10)
dChip.clustering(data.quantiles[vargenes_0.5_10,],annotation,main="",col_size=1,p_classes=samplesheet$"Replicate")
classes <- samplesheet$"Replicate"
perform.pca(data.quantiles[vargenes_0.5_10,],classes,html=T,outputfile="pca_vargenes_0.5_10")
# legend positions "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center"
perform.pca(data.quantiles[vargenes_0.5_10,],classes,html=F,outputfile="pca_vargenes_0.5_10",pos_legend="bottomright")
perform.pca(data.quantiles[vargenes_0.5_10,],classes,html=F,outputfile="pca_vargenes_0.5_10",plot_3d=F)

# DE genes
background <- getBackground(data.quantiles,raw.calls)
print(paste("Background signal intensity:",background))
DEgenes_A_vs_B <- getDEgenes(data.frame(data.quantiles[,which(classes=="treatmentA")]),
        data.frame(data.quantiles[,which(classes=="treatmentB")]),annotation,outputfile="DEgenes_A_vs_B",fc=2,diff=background,pval=0.05)
length(DEgenes_A_vs_B)
dChip.clustering(data.quantiles[DEgenes_A_vs_B,],annotation,main="",col_size=1,p_classes=samplesheet$"Replicate")
perform.pca(data.quantiles[DEgenes_A_vs_B,],classes,html=T,outputfile="pca_DEgenes")

# ANOVA => compares >2 classes
ANOVAgenes <- getANOVAgenes(data.quantiles,classes,fc=2,diff=background,pval = 0.05,annotation,outputfile="ANOVAgenes")
length(ANOVAgenes)
dChip.clustering(data.quantiles[ANOVAgenes,],annotation,main="",col_size=1,p_classes=samplesheet$"Replicate")
perform.pca(data.quantiles[ANOVAgenes,],classes,html=T,outputfile="pca_ANOVAgenes")

# save complete data => change annotation
colnames(data.quantiles) <- samplesheet$"Group"
writetable(cbind(annotation[rownames(data.quantiles),],data.quantiles),"data_all")


