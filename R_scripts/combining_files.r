fn=function(file)
{
file_list <- unlist(strsplit(file,split=","))
setwd("D:/combine")
for (file in file_list)
{
dataset=read.delim(file, header=TRUE, sep="\t")
dataset=cbind(dataset,dataset)                     
}
names(dataset)[1]=paste("Probeid")
dataset=dataset[,-match(c("ProbeID"),names(dataset))]

                                                                   
write.table(dataset,file="D:/output/illumina.txt",sep="\t",col.names=NA)

return ("illumina.txt")
}
