library(riboSeqR)
library(Rsamtools)
fastaCDS=findCDS("~/Transcriptomes/smed_dd_v6.fasta",startCodon = c("ATG"),stopCodon = c("TAG","TAA","TGA"))
riboseqR_ORFs=as.data.frame(fastaCDS)
rpf_counts=read.delim("/home/tsekara/Novel_ORFs_Prediction/time_series/vs_dd_v6_transcriptome/RPF_0hrs/linker_free_RPFs/Pooled_sam/sam2profile_Pooled_0hrs_RPF.txt")
translated_riboseq_ORFs=riboseqR_ORFs[riboseqR_ORFs$seqnames %in% rpf_counts$Transcript_Id,]
write.table(translated_riboseq_ORFs,"~/Translated_transcripts_riboseqr_orfs/Translated_trnascripts_riboseeqR_orfs.txt",sep="\t",quote=F,row.names=F)

########################################### in Linux ########################################
mkdir Translated_transcripts_riboseqr_orfs
cd Translated_transcripts_riboseqr_orfs
cp Translated_transcripts_riboseqr_orfs.txt ~/Translated_transcripts_riboseqr_orfs
awk '{ print >> $1".txt" }' Translated_transcripts_RiboseqR_ORFS.txt

########################################## in R ##############################################
setwd("~/Translated_transcripts_riboseqr_orfs/")
############## filter trasncripts with raw counts >=10 ###############################
gfp_1=rpf_counts[rpf_counts$gfp_1>=10,]
################### Delete the rRNA ##############################
gfp_1=gfp_1[-c(1,2),]
################### Extractingt he Transcript_Ids only ##############################
length(gfp_1$Transcript_Id)
pooled_0hrs_ids=paste(gfp_1$Transcript_Id,".txt",sep="")
head(pooled_0hrs_ids)
################### Read the input Bam file #############################
riboDat_gfp_1=readRibodata("~/path/to/Sorted_file.bam",replicate="GFP_1")

################### load the libraries for utilising multi-cores and registering them ############################
library("parallel")
library("foreach")
library("doParallel")

cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl, cores = detectCores() - 1)

########################## Run the for each loop  for itertaing through every txt file in folder for counting the frames and plotiing ############################
foreach(i=1:length(pooled_0hrs_ids),.packages = c("riboSeqR"),.combine=rbind) %dopar%
{
  cds=read.delim(pooled_0hrs_ids[i],as.is=T,header=F)
  
  cds=GRanges(cds$V1,IRanges(start=cds$V2,end=cds$V3))
  cds$frame=(start(cds)-1)%%3
  fCs=frameCounting(riboDat_pooled_0hrs,cds,lengths = 27:32)
  fS=readingFrame(rC=fCs,lengths = 27:32)
  write.table(fS,file=paste("/home/tsekara/Novel_ORFs_Prediction/time_series/vs_dd_v6_transcriptome/RPF_0hrs/linker_free_RPFs/Pooled_sam/Codon_periodicity/Pooled_0hrs_frame_counting_matrix/",pooled_0hrs_ids[i],sep=""),sep="\t",quote=F,row.names=F)
  pdf(file=paste("/home/tsekara/Novel_ORFs_Prediction/time_series/vs_dd_v6_transcriptome/RPF_0hrs/linker_free_RPFs/Pooled_sam/Pooled_0hrs_Plots/",gsub("txt","pdf",pooled_0hrs_ids[i]),sep=""))
  plotFS(fS)
  dev.off()
}