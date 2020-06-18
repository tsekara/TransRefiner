###### Thileepan Sekaran Pipeline to compute DEGs using limma and voom ############
                                            
################### L & Voom ###################
######### Loading reauired libraries ##############
library("limma")
library("edgeR")
##### Assuming the counts generated samtools idxstats. The count data(idxcounts.txt) file for each sample can be integrated into a matrix file using the following script ########  
counts.df=NULL
experiment.names <- c("control_1","control_2","control_3","treated_A_1","treated_A_2","treated_A_3","treated_B_1","treated_B_2","treated_B_3");
for(expName in experiment.names)
{
  reads.df <- read.delim(sprintf("%s_idxcounts.txt",expName),col.names = c("refName","refLength","mapped","unmapped"));
  counts.df <- cbind(counts.df, reads.df$mapped);
}
rownames(counts.df) <- reads.df$refName;
colnames(counts.df) <- experiment.names;
########### samplesheet creation ####################    
samplesheet=data.frame(row.names=colnames(counts.df),condition=c("control","control","control","treated_A","treated_A","treated_A","treated_B","treated_B","treated_B"))
######## Creating Desing matrix ########################    
design=model.matrix(~0+samplesheet$condition)
colnames(design)=c("control","treated_A","treated_B")
###################### Creating a Contrast matrix #######################
cont.matrix=makeContrasts("control-treated_A","control-treated_B",levels=design)
############ Filtering very low counts whose sum of counts is less than 10 ################      
isexp=rowSums(cpm(counts.df) > 10) >= 2
counts.df <- counts.df[isexp,]
############ Creating the DGE list from count data #####################    
dge=DGEList(counts=counts.df)
dge <- calcNormFactors(dge)
v <- voom(dge,design,plot=TRUE)
######## LIMMA pipeline for DEGs for only two conditions ######################
fit=lmFit(v,design)
fit=eBayes(fit)    
topTable(fit, adjust="BH")
############ LIMMA pipeline for DEGs for 2 or more conditions ########################
fit=lmFit(v,design)
fit2=contrasts.fit(fit,cont.matrix)
fit2=eBayes(fit2)    
topTable(fit2,coef="control-treated_A",adjust="BH") ####### Toptable for condition1 #######
topTable(fit2,coef="control-treated_B",adjust="BH") ####### Toptable for condition2 #######    
toptable(fit2,adjust="BH") #################### Toptable for both condition ###########

results= decideTests(fit2) ################### It assigns the outcome of every condition #####  
vennDiagram(results)

######### to find ANOVA and display top 30 genes ###################
topTableF(fit2, number=30)

  
    
    
    
    
    
    
    
    
    
    
    
    
    
