setwd("~/XXX")
library( DESeq2 )
rpf_deg<-read.table("xx.txt", fill=T,header=T, sep="\t")
head(RPF)
mrna_deg<-read.table("yy.txt", fill=T,header=T, sep="\t")
head(MRNA)
RPF_MRNA<-merge(rpf_deg, mrna_deg, by.x="idR", by.y="idM", all=T)
write.table(RPF_MRNA, "xyxy.csv",sep=",", col.names=NA, row.names=T )
data.df<-read.table("rpf_mrna.csv", fill=T,header=T, row.names=1, sep=",")

head(data.df)
data.sub <- subset(data.df, data.df$baseMeanR>150)
head (data.sub)
data.sub2 <- subset(data.sub, data.sub$baseMeanM>150)
head (data.sub2)

summary(lm(data.sub2$log2FoldChangeR~data.sub2$log2FoldChangeM))

dim(data.sub2[ data.sub2$padjM < .01 & data.sub2$padjR< .01, ])
data.sub2Sig01 <- data.sub2[ data.sub2$padjM < .01 & data.sub2$padjR< .01, ]

dim(data.sub2[ data.sub2$padjM < .05 & data.sub2$padjR< .05, ])
data.sub2Sig05 <- data.sub2[ data.sub2$padjM < .05 & data.sub2$padjR< .05, ] 

dim(data.sub2[ data.sub2$pvalR < .1 & data.sub2$pvalM< .1, ])
data.sub2Sig1 <- data.sub2[ data.sub2$pvalR < .1 & data.sub2$pvalM< .1, ]

write.csv(data.sub2Sig01, "xyxy_sig.csv")

library ("scales")




mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0,1, 0, 0))
plot(data.sub2$log2FoldChangeR~data.sub2$log2FoldChangeM,  family="sans", font.axis=2, cex.axis=1.4, cex.lab=1.6, font.lab=2, xlab="log2 foldchange mRNA", ylab="log2 foldchange RPF", pch=19, cex=.25, xlim=c(-6, 6), ylim=c(-6,6), col = ifelse(data.sub2$padjM< .1 & data.sub2$padjR< .1, "red", ifelse(data.sub2$padjM< .1, "dark green", ifelse(data.sub2$padjR< .1, "orange", "grey" ) ) ))
head(data.sub2)
RPF_MRNA<-merge(res_rpf, res_mrna, by.x="idR", by.y="idM", all=T)
data.df=RPF_MRNA
data.sub <- subset(data.df, data.df$baseMeanR>50)
head(data.sub)
data.sub2 <- subset(data.sub, data.sub$baseMeanM>50)
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0,1, 0, 0))
plot(data.sub2$log2FoldChangeR~data.sub2$log2FoldChangeM,  family="sans", font.axis=2, cex.axis=1.4, cex.lab=1.6, font.lab=2, xlab="log2 foldchange mRNA", ylab="log2 foldchange RPF", pch=19, cex=.25, xlim=c(-6, 6), ylim=c(-6,6), col = ifelse(data.sub2$padjM< .1 & data.sub2$padjR< .1, "red3", ifelse(data.sub2$padjM< .1, "coral", ifelse(data.sub2$padjR< .1, "blue", "gray50" ) ) ))

selectedr.genes<-ifelse(data.sub2$padjM< .01 & data.sub2$padjR< .01, FALSE, ifelse(data.sub2$padjM< .01, FALSE, ifelse(data.sub2$padjR< .01, TRUE, FALSE ) ) )
selectedr.genes[is.na(selectedr.genes)]<-FALSE
dim(data.sub2[ selectedr.genes,])
only.RPF<-data.sub2[ selectedr.genes,]
write.csv(only.RPF, "xyxy_RPF.csv")

selectedm.genes<-ifelse(data.sub2$padjM< .01 & data.sub2$padjR< .01, FALSE, ifelse(data.sub2$padjM< .01, TRUE, ifelse(data.sub2$padjR< .01, FALSE, FALSE ) ) )
selectedm.genes[is.na(selectedm.genes)]<-FALSE
dim(data.sub2[ selectedm.genes,])
only.mRNA<-data.sub2[ selectedm.genes,]
write.csv(only.mRNA, "xyxy_mRNA.csv")

selectedmr.genes<-ifelse(data.sub2$padjM< .01 & data.sub2$padjR< .01, TRUE, ifelse(data.sub2$padjM< .01, FALSE, ifelse(data.sub2$padjR< .01, FALSE, FALSE ) ) )
selectedmr.genes[is.na(selectedmr.genes)]<-FALSE
dim(data.sub2[ selectedmr.genes,])
only.overlap<-data.sub2[ selectedmr.genes,]
write.csv(only.overlap, "xyxy_mrna_rpf.csv")