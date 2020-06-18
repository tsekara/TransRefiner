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



      raw.Data=read.delim(file="D:/output/illumina.txt",row.names = 1, dec = ".")
      raw.expression <- raw.Data[,seq(1,dim(raw.Data)[2],2)]
      dim(raw.expression)   
      raw.calls <- raw.Data[,seq(2,dim(raw.Data)[2],2)]
      dim(raw.calls)

      library (limma)
      data.quantiles<- normalizeBetweenArrays (as.matrix (raw.expression), method = "quantile")
      colnames (data.quantiles) <- colnames (raw.expression)
      rownames (data.quantiles) <- rownames (raw.expression)   
      


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
       
    #writetable(variable.genes,"Variable_genes",p_row=T)
    
    write.table(variable.genes,file ="D:/output/Variable_genes.txt",sep="\t",col.names=T,row.names=T)
    
    var=sum(variable.genes)
  
  
    return("Variable_genes.txt")
}
