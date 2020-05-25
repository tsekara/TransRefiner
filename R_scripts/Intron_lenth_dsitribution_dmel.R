intron_len=read.delim("Dmel_intronseq.fasta.fai",header=F)
head(intron_len)
intron_len$V3=NULL
intron_len$V4=NULL
intron_len$V5=NULL
plot(log(intron_len))
rownames(intron_len)=intron_len$V1
intron_len$V1=NULL
head(intron_len)
plot(log2(intron_len$V2))
plot(density(log2(intron_len$V2)))
log2(70)
plot(density(log2(intron_len$V2)),main = "Drosophila_Intron_length_Distribution",xlab = "Log2_intron_len")
