combine=function(fnames)
{
strp=unlist(strsplit(fnames,split=","))
n= numeric(length(strp))
setwd("D:/combine")
for(i in seq_along(strp))
{
    d = read.delim(strp[i])
    n[i] = dim(d)[1]
}

setwd("D:/output")
pdf(file=paste("Probe_Ids_comparision.pdf",sep=""))
plot(n,ylab=" Number of Probe IDs",col=rainbow(n))
dev.off()
return("Probe_Ids_comparision.pdf")
}
