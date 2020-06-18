
            combine=function(file)
            {
    
    
      split_list <- unlist(strsplit(file,split=","))
    	setwd("E:/combine")
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
    
       samplesheet <- read.delim(""H:/illumina_test/SampleSheet.txt")
      
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
  
      
              ##################Variable Genes###################

  getsdvsmean <- function(vec)

  {
    sdvsmean <- sd(vec,na.rm=T)/mean(vec,na.rm=T)
    return(sdvsmean)
  }

  is.variable <- function(quot)
  {
    if (0.5 < quot && quot < 10)
    {
      return(TRUE)
    }
    else
    {
      return(FALSE)
    }
  }

  data.mean <- apply(data.quantiles,1,mean,na.rm=T)
  data.std <- apply(data.quantiles,1,sd,na.rm=T)
  quot <-data.std/data.mean

  variable.genes <- sapply(quot,is.variable)
  var.genes=sum(variable.genes)


     classes = colnames(data.quantiles[variable.genes,])

                ################dchipClustering#######################
require (amap)
require(marray)
zscore<-function (row)
{
mw<-mean (row, na.rm= TRUE)
std<-sd (row, na.rm= TRUE)
row<- (row-mw)/std
return (row)
}
col.pal<-maPalette (low = "darkblue", high = "darkred", mid = "white", k = 1000)
############################# standardize rows ################################
p_zscore = FALSE
if(p_zscore==TRUE)
{
clusterdata<- t(apply (data.quantiles[variable.genes,], 1, zscore))
rownames (clusterdata) <- annotation[as.character(rownames (data.quantiles[variable.genes,])), 1]
}
else
{
clusterdata<- as.matrix(data.quantiles[variable.genes,])
rownames(clusterdata) <- annotation[as.character(rownames(data.quantiles[variable.genes,])),1]
}
########################### cluster data #########################################
clustersamples<-as.dendrogram(hcluster (t (clusterdata), method = "correlation", link ="average"))
clustergenes<-as.dendrogram(hcluster (clusterdata, method ="correlation" , link = "average"))
col.pal<-maPalette (low = "darkblue", high = "darkred", mid = "white", k = 1000)
csc<-rainbow (max (max (as.integer(factor (colnames(data.quantiles[variable.genes,]))))))[as.integer( factor (colnames(data.quantiles[variable.genes,])))]
setwd("D:/output")
pdf(file=paste("dchipclustering_var_genes.pdf",sep=""))
heatmap (clusterdata, Rowv = clustergenes, Colv = clustersamples, col = col.pal, main = "Clustering", labRow = rownames (clusterdata), labCol = classes, cexCol = 0.3, cexRow = 0.3, ColSideColors = csc, margins = c(6,6) )   
detach("package:amap")
 dev.off()
 return("dchipclustering_var_genes.pdf")
}


