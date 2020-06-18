
    combine=function(file)
    {
        split_list <- unlist(strsplit(file,split=","))
      	setwd("D:/combine")
	     dataset <- do.call("cbind",lapply(split_list,FUN=function(files)
	         {
              read.table(files,header=TRUE, sep="\t")
	         }
	                                                           )
	                                                           )
	    names(dataset)[1]=paste("Probe_ID")
      drop=c("ProbeID")
      dataset=dataset[,!(names(dataset)%in%drop)]

	     write.table(dataset,file ="D:/output/illumina.txt",sep="\t",col.names=T,row.names=F)

      
    
              #########  reading illumina file   ###########
   raw.Data <- read.delim(file ="D:/output/illumina.txt", row.names = 1, dec = ".") 
   raw.expression <- raw.Data[,seq(1,dim(raw.Data)[2],2)]
   raw.calls <- raw.Data[,seq(2,dim(raw.Data)[2],2)]
  
  
          ###################  Sample Sheet ####################
  
   samplesheet <- read.delim("H:/illumina_test/SampleSheet.txt")
   
   colnames(raw.expression) <- samplesheet$"Group"
   colnames(raw.calls) <- samplesheet$"Group"

         #########  Normalisation  ############
   
    library (limma)
    data.quantiles<- normalizeBetweenArrays (as.matrix (raw.expression), method = "quantile")
    colnames (data.quantiles) <- colnames (raw.expression)
    rownames (data.quantiles) <- rownames (raw.expression)
                                                               
        #########  Reading Annotation file  ############
   annotation=read.delim (file="H:/illumina_test/MouseWG-6_V2.txt", row.names = "Array_Address_Id", dec = ",")
   
   classes <- samplesheet$"Replicate"
   sample1=data.frame(data.quantiles[,which(classes=="treatmentA")])
   sample2=data.frame(data.quantiles[,which(classes=="treatmentB")])
   
        #########  Calculating Differential Expressed genes   ############
   
   fc=2
   pval=0.05
   diff=100
   output_all=F

  library(multtest)
  library(stats)
  mean_s1 <- apply(sample1,1,mean)
  mean_s2 <- apply(sample2,1,mean)
  sd_1 <- apply(sample1,1,sd)
  sd_2 <- apply(sample2,1,sd)
  
  ######################## isfoldchange #################################
  
 isfoldchange= function(sample,fc)
{
  change.fc <- function(value)
  {
    if (value < 1)
    {
      new_value <- -(1/value)
    }
    else
    {
      new_value <- value
    }
    return(new_value)
  }
  
  isdifferent <- 0
  fold <- sample[1]/sample[2]
  if (is.na(fold))
  {
    fold<-1
  }
  if (fold > fc || fold < (1/fc))
  {
    isdifferent <- 1
  }
  result <- c(sapply(fold,change.fc),isdifferent)
  return(result)
}

                                            
  ######################## Getfoldchange ################################
  

  
   getfoldchange = function(sample1,sample2,fc)
    {
          t(apply(as.matrix(cbind(sample1,sample2)),1,isfoldchange,fc))
    }   
  
  fc_s1_vs_s2 <- getfoldchange(mean_s1,mean_s2,fc)
  
  teststat_s1_vs_s2 <- mt.teststat(cbind(sample1,sample2),c(rep(0,dim(sample1)[2]),rep(1,dim(sample2)[2])))
  pvals_s1_vs_s2 <- 2*(1-pnorm(abs(teststat_s1_vs_s2)))
  pval_adjust <- p.adjust(pvals_s1_vs_s2, method = "BH")
  difference <- abs(mean_s1-mean_s2)

  DEgenes_index <- intersect(intersect(which(fc_s1_vs_s2[,2]==1),which(pval_adjust<pval)),which(difference>diff))

  if (output_all == T)
  {
    is_diff <- rep(0,length(rownames(sample1)))
    is_diff[DEgenes_index] <- 1
    output <- cbind(annotation[rownames(sample1),],mean_s1=round(mean_s1,2),mean_s2=round(mean_s2,2),sd_1=sd_1,sd_2=sd_2,FC=round(fc_s1_vs_s2[,1],2),pvalue=round(pval_adjust,4),diff=round(difference,2),is_DE=is_diff)
  }
  else
  {
    output <- cbind(annotation[rownames(sample1[DEgenes_index,]),],mean_s1=round(mean_s1[DEgenes_index],2),mean_s2=round(mean_s2[DEgenes_index],2),sd_1=sd_1[DEgenes_index],sd_2=sd_2[DEgenes_index],FC=round(fc_s1_vs_s2[DEgenes_index,1],2),pvalue=round(pval_adjust[DEgenes_index],4),diff=round(difference[DEgenes_index],2))
  }
    setwd("D:/output")
  
    #writetable(output,outputfile,p_row=T)
  
    write.table(output,file ="D:/output/DE_genes.txt",sep="\t",col.names=T,row.names=T)
  
    return("DE_genes.txt")
}
