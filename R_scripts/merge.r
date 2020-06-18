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

return ("illumina.txt")
}