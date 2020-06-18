        try=function(file)
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

        ######Reading the merged illumina text files along with annotation file and SampleSheet#############
    
       samplesheet <- read.delim("H:/illumina_test/SampleSheet.txt")
      
       classes <- samplesheet$"Replicate"
      
       annotation=read.delim (file="H:/illumina_test/MouseWG-6_V2.txt", row.names = "Array_Address_Id", dec = ",")

      raw.Data=read.delim(file="D:/output/illumina.txt",row.names = 1, dec = ".")
      raw.expression <- raw.Data[,seq(1,dim(raw.Data)[2],2)]
      dim(raw.expression)   
      raw.calls <- raw.Data[,seq(2,dim(raw.Data)[2],2)]
      dim(raw.calls)

      ##########################Normalization##################################

      library (limma)
      data.quantiles<- normalizeBetweenArrays (as.matrix (raw.expression), method = "quantile")
      colnames (data.quantiles) <- colnames (raw.expression)
      rownames (data.quantiles) <- rownames (raw.expression)   
      
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
    
    
    getBackground=function(norm_data,calls)
      {
          intensity_signal <- vector()
          detection_pval <- vector()
          for(i in 1:ncol(calls))
      {
          intensity_signal <- c(intensity_signal,norm_data[,i])
          detection_pval <- c(detection_pval,calls[,i])
      }
          temp <- intersect(which(detection_pval<0.055),which(detection_pval>0.045))
          background_model <- lm(intensity_signal[temp]~detection_pval[temp])
          background <- background_model$coefficients[1]+background_model$coefficients[2]*0.05
          return(background)
      }
      background <- getBackground(data.quantiles,raw.calls)
      
      
 getDEgenes=function(sample1,sample2,annotation,fc=2,pval=0.05,diff=100,output_all=F,outputfile="DEgenes")
{
  library(multtest)
  library(stats)
  mean_s1 <- apply(sample1,1,mean)
  mean_s2 <- apply(sample2,1,mean)
  sd_1 <- apply(sample1,1,sd)
  sd_2 <- apply(sample2,1,sd)
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
  writetable(output,outputfile,p_row=T)

  return(rownames(sample1[DEgenes_index,]))
  return(length(DEgenes_index))
}
    DEgenes_A_vs_B <- getDEgenes(data.frame(data.quantiles[,which(classes=="treatmentA")]),
        data.frame(data.quantiles[,which(classes=="treatmentB")]),annotation,outputfile="DEgenes_A_vs_B",fc=2,diff=background,pval=0.05)
        
        
dchip.clustering<-function(data,annotation,p_method="correlation",p_link="average",p_zscore=F,p_heatmap=T,main_annot="Clustering",classes=colnames(data),
col_size=0.3,row_size=0.3,genes=rownames(data),p_classes = colnames(data), c_margins = c(6,6))
{
  require(amap)
  require(marray)

  zscore <- function(row)
  {
    mw <- mean(row,na.rm=T)
    std <- sd(row,na.rm=T)
    row <- (row-mw)/std
    return(row)
  }
  col.pal <- maPalette(low="darkblue", high="darkred",mid="white",k=1000)

  if (p_zscore==TRUE)
  {
    clusterdata <- t(apply(data,1,zscore))
    rownames(clusterdata) <- annotation[as.character(rownames(data)),1]
  }
  else
  {
    clusterdata <- as.matrix(data)
    rownames(clusterdata) <- annotation[as.character(rownames(data)),1]
  }
    clustersamples <- as.dendrogram(hcluster(t(clusterdata),method=p_method,link=p_link))
    clustergenes <- as.dendrogram(hcluster(clusterdata,method=p_method,link=p_link))
    col.pal <- maPalette(low="darkblue", high="darkred",mid="white",k=1000)
    csc <- rainbow(max(max(as.integer(factor(p_classes)))))[as.integer(factor(p_classes))]
    heatmap(clusterdata,Rowv=clustergenes,Colv=clustersamples,col=col.pal,main=main_annot,labRow=rownames(clusterdata),labCol=classes,cexCol=col_size,cexRow=row_size,ColSideColors=csc, margins = c_margins)
    dev.off()
  detach("package:amap")
}
 
 dChip.clustering(data.quantiles[DEgenes_A_vs_B,],annotation,main="",col_size=1,p_classes=samplesheet$"Replicate")
 
 return("dchip.clustering_DEgenes.pdf")
 }