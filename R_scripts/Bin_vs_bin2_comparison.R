
cds_regions=read.delim("extended_5pp_3531seqs_CDS.bed",header=F)
cds_regions$V5=NULL
Bin1_regions=cds_regions %>% mutate(Bin1_Start=ifelse((V2-13 >0 & V4=="+" ),V2-13,V3-13)) %>% mutate(Bin1_End=ifelse((V2-13 > 0 & V4=="+"),V2+37,V3+37))
write.table(Bin1_regions,"Bin1_CDS_-13nts_37nts.bed",sep="\t",quote=F,col.names = F,row.names = F)

Bin2_regions=cds_regions %>% mutate(Bin1_Start=ifelse((V2-13 >0 & V4=="+" ),V2-13,V3-13)) %>% mutate(Bin1_End=ifelse((V2-13 > 0 & V4=="+"),V2+37,V3+37)) %>% mutate(Bin2_start=ifelse(V4=="+",Bin1_End+1,Bin1_Start-51)) %>%   mutate(Bin2_end=ifelse((V4=="+"),Bin2_start+50,Bin2_start+50)) %>% select(V1,Bin2_start, Bin2_end)
write.table(Bin2_regions,"Bin2_CDS_50nts.bed",sep="\t",quote=F,col.names = F,row.names = F)


bin1_new=read.delim("Bin1_counts.txt",header=F)
colnames(bin1_new)[4]="Bin1_counts"
bin2_new=read.delim("Bin2_counts.txt",header=F)
colnames(bin2_new)[4]="Bin2_counts"
bin1_bin2_new=merge(bin1_new,bin2_new,by="V1")
bin1_bin2_new %>% filter(Bin1_counts>20 & Bin2_counts >20) %>% mutate(Ratio=Bin1_counts/Bin2_counts) %>% filter(Ratio>2) %>% dim()
