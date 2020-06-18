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



      raw.Data=read.delim(file="D:/output/illumina.txt",row.names = 1, dec = ".")
      raw.expression <- raw.Data[,seq(1,dim(raw.Data)[2],2)]
      dim(raw.expression)
      raw.calls <- raw.Data[,seq(2,dim(raw.Data)[2],2)]
      dim(raw.calls)

  setwd("D:/output")
  directory=getwd()

  IDs <- colnames(raw.expression)


 pdf(file="Unnormalised_Scatterplots.pdf",onefile=TRUE)
  for (i in 1:(dim(raw.expression)[2]-1))
  {
    for( j in i:(dim(raw.expression)[2]) )
    {

      if (i != j)
      {

        correlation <- round(cor(raw.expression[,i],raw.expression[,j]),2)
        maximum <- max(log2(raw.expression[,i]))
        minimum <- min(log2(raw.expression[,i]))
        plot(log2(raw.expression[,i]),log2(raw.expression[,j]),xlab=IDs[i],ylab=IDs[j],pch='.',text(maximum-2,minimum+0.5,labels=paste("R = ",correlation,sep=""),pos=4,offset=0))

      }

    }

  }

   dev.off()
   return("Unnormalised_Scatterplots.pdf")
}