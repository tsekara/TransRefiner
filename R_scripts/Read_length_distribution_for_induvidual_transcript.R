# R and Linux script to compute read distribution based on its length per transcript.
# The linux code must be run on the terminal and the output text file must be loaded into R



################# Linux awk code #####################################################
samtools view bam_file | gawk '(! and(4, $2))' | awk '{print $3, length($10)}' > read_length_distribution.txt




################################# R-code ###################################################
test=read.delim("read_length_distribution.txt",header = F)
names(test)=c("ID","Read_Length")
head(setDT(test)[,list(len_26=sum(Read_Length==26), len_27=sum(Read_Length==27),len_28=sum(Read_Length==28),len_29=sum(Read_Length==29),len_30=sum(Read_Length==30),len_31=sum(Read_Length==31),len_32=sum(Read_Length==32),len_33=sum(Read_Length==33),len_34=sum(Read_Length==34)),by="ID"])
read_length_distribution_per_transcript=setDT(test)[,list(len_26=sum(Read_Length==26), len_27=sum(Read_Length==27),len_28=sum(Read_Length==28),len_29=sum(Read_Length==29),len_30=sum(Read_Length==30),len_31=sum(Read_Length==31),len_32=sum(Read_Length==32),len_33=sum(Read_Length==33),len_34=sum(Read_Length==34)),by="ID"]
write.table(read_length_distribution_per_transcript,"Read_length_distribution.txt",sep="\t",quote=F,row.names=F)
