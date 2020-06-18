

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

	write.table(dataset,file ="E:/output/illumina.txt",sep="\t",col.names=T,row.names=F)



    raw.Data=read.delim(file="E:/output/illumina.txt",row.names = 1, dec = ".")
    raw.expression <- raw.Data[,seq(1,dim(raw.Data)[2],2)]
    dim(raw.expression)
    raw.calls <- raw.Data[,seq(2,dim(raw.Data)[2],2)]
    dim(raw.calls)

      library (limma)
      data.quantiles<- normalizeBetweenArrays (as.matrix (raw.expression), method = "quantile")
      colnames (data.quantiles) <- colnames (raw.expression)
      rownames (data.quantiles) <- rownames (raw.expression)   
      

    samplesheet=read.delim("E:/output/ss.txt")

    setwd("E:/output")
    present <- vector()
    for(i in 1:ncol(raw.calls))
    {
      present[i] <- sum(raw.calls[,i]<0.05)/nrow(raw.calls)*100
    }
    names(present)= samplesheet[,1]
    write.table(present,"Present_call.txt",row.names=T,col.names="ArrayIDs",sep="\t")
    return("Present_call.txt")
  }