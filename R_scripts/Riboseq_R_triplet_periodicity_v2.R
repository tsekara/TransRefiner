################################# Codon periodicity with RiboseqR ###############################

# Though RiboseqR has inbuilt function to view the triplet # periodicity but peridicity only for a particular 
# length can visulaized at a time. Therefore I modified the riboSeqR sligthly so that triplet peridicity can 
# be viewed for multiple length(27:32) at a time can be visualized as box plot as well export the frame
# counts in text file which can be used for further down stream analysis.


library(riboSeqR)
library(GenomicAlignments)

#################### Load the RPF-count table inside R ########################################### 
rpf_counts=read.delim("/home/tsekara/Novel_ORFs_Prediction/time_series/vs_dd_v6_transcriptome/RPF_0hrs/linker_free_RPFs/Pooled_sam/sam2profile_Pooled_0hrs_RPF.txt")
#################### Retain transcripts with RPF-counts > X(10 or 20 or 30) #####################################
rpf_counts=rpf_counts[rpf_counts$Reads >10,]
#################### remove the rRNA counts ######################################################
rpf_counts=rpf_counts[-c(1,2),]
########################### Get the transcripts IDs #############################################
ids=rpf_counts$Transcript_Id
########################### adding suffix ".txt" to access the all_frames_ORFs text files ##################
files=paste(ids,".txt",sep="")
########################## setting "riboSeqR_all_frames_ORFs" directory as working directory ########################
setwd("~/riboseqR_all_frame_CDS/")

riboseq_R_triplet_periodicty=function(ids,files,frame_count_path,frame_plot_path)
{  
for(i in 1:length(ids))
  {
  #cds=read.delim(files[i],as.is=T,header=F)
  cds=read.delim(paste("~/riboseqR_all_frame_CDS/",files[i],sep=""),as.is = T,header = F)
  cds=GRanges(cds$V1,IRanges(start=cds$V2,end=cds$V3))
  cds$frame=(start(cds)-1)%%3
  system(paste0('samtools view -F 4 -@ 40 -bS -o ~/output.bam ~/Transcriptome_Enhancement/Pooled/Sorted_pooled_master.bam',sep=" ",sub_ids[i]))
  riboDat=readRibodata("~/output.bam",replicates="Pooled")
  fCs=frameCounting(riboDat,cds,lengths = 27:32)
  fS=readingFrame(rC=fCs,lengths = 27:32)
  write.table(fS,file=paste("/home/tsekara/",ids[i],".txt",sep=""),sep="\t",quote=F,row.names=F)
  #pdf(file=paste("/home/tsekara/",ids[i],".pdf",sep=""))
  pdf(file=paste(frame_plot_path,ids[i],".pdf",sep=""))
  plotFS(fS)
  dev.off()
  }
}  