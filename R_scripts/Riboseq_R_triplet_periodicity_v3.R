library(riboSeqR)
library(GenomicAlignments)
library(doParallel)

##################### Loading a Sorted Bam file ##########################
riboDat=readRibodata("~/path/to/Sorted_file.bam",replicate="Pooled")

##################### Load the FASTA file #################################

fastaCDS=findCDS("~/Path_to/transcripts.fasta",startCodon = c("ATG"),stopCodon = c("TAG","TAA","TGA"))

##################### SPLIT the above Fasta on IDs as list #################################################

riboseqR_ORFs=as.data.frame(fastaCDS)
induvidual_transcripts <- split( riboseqR_ORFs , f = riboseqR_ORFs$seqnames )

##############################################  Method 1 ####################################################

if((dir.exists("./Frame_Plots/")==FALSE)&&dir.exists("./Frame_Counts/")==FALSE)
{
  dir.create(path = "./Frame_Plots/")
  dir.create(path = "./Frame_Counts/")
}

foreach(i=1:length(induvidual_transcripts),.packages = c("riboSeqR"),.combine=rbind) %dopar%
{
  cds=as.data.frame(induvidual_transcripts[i])
  cds=GRanges(cds[,1],IRanges(start=cds[,2],end=cds[,3]))
  cds$frame=(start(cds)-1)%%3
  fCs=frameCounting(riboDat,cds,lengths = 27:32)
  fS=readingFrame(rC=fCs,lengths = 27:32)
  write.table(fS,file=paste("./Frame_Counts/",names(induvidual_transcripts[i]),".txt",sep=""),sep="\t",quote=F,row.names=F)
  pdf(file=paste("./Frame_Plots/",names(induvidual_transcripts[i]),".pdf",sep=""))
  plotFS(fS)
  dev.off()
}

###########################################  Method 2 #############################################################

if((dir.exists("./Frame_Plots/")==FALSE) && dir.exists("./Frame_Counts/")==FALSE)
{
  dir.create(path = "./Frame_Plots/")
  dir.create(path = "./Frame_Counts/")
}

  for(i in 1:length(induvidual_transcripts))
  {
    cds=as.data.frame(induvidual_transcripts[i])
  # cds=read.delim(paste("~/riboseqR_all_frame_CDS/",files[i],sep=""),as.is = T,header = F)
    cds=GRanges(cds[,1],IRanges(start=cds[,2],end=cds[,3]))
    cds$frame=(start(cds)-1)%%3
    system(paste0('samtools view -F 4 -@ 40 -bS -o ~/output.bam /home/tsekara/Transcriptome_Enhancement/Transdecoder_test/extended_transcripts/Sorted_extended_transcripts.bam',sep=" ",names(induvidual_transcripts)[i]))
    riboDat=readRibodata("~/output.bam",replicates="Pooled")
    fCs=frameCounting(riboDat,cds,lengths = 27:32)
    fS=readingFrame(rC=fCs,lengths = 27:32)
    write.table(fS,file=paste("./Frame_Counts/",names(induvidual_transcripts[i]),".txt",sep=""),sep="\t",quote=F,row.names=F)
    #pdf(file=paste("/home/tsekara/",induvidual_transcripts[i],".pdf",sep=""))
    pdf(file=paste("./Frame_Plots/",names(induvidual_transcripts[i]),".pdf",sep=""))
    plotFS(fS)
    dev.off()
  }