


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
  

      


       
      samplesheet=read.delim("H:/illumina_test/SampleSheet.txt")
      classes <- samplesheet$"Replicate"
      annotation=read.delim ( "H:/illumina_test/MouseWG-6_V2.txt", row.names = "Array_Address_Id", dec = ",")

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
  var.genes=sum(variable.genes)
  
    require(pcurve)
    require(made4)
    dataset_pca<- pca (t (data.quantiles[variable.genes,]))
    pcadata<- as.matrix(cbind (dataset_pca$pcs[,1], dataset_pca$pcs[,2], dataset_pca$pcs[,3]))
    setwd("E:/output")
        if (max (abs (pcadata)) > 100000)
            {
                pcadata<-pcadata/2
            }
        if (TRUE)
            {
                html3D (pcadata, classvec = classes, writehtml = TRUE, filenamebase = "var_genes_PCA")
            }
          else
            {                                   
                do3d (pcadata, angle = 40, classvec = classes, cex.symbols = 1)               
            }
                     
                  return("var_genes_PCA.html")
}