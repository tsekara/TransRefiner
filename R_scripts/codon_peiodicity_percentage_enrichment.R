################## Validation of Extended transcripts with frame Enrichment ##################

colMax <- function(data) sapply(data, max, na.rm = TRUE) # function to find the max value in columns

# Copy the folder "Frame_counting matrix" in the same directory and rename it as "Frame_counting_without_frame_info" and copy the script "removing_last_line_of_every_text_file.sh"and
# execute it so that the frame information(last line) in every text file is removed.

# Load the file names of frame_counting matrix generated from custom_riboseqR function

file_names=dir(pattern = "*.txt")

# Function for retaining only the maximum value 

frame_percentage=function(file_names,temp_path)
{  
for (i in 1:length(file_names))
#foreach(i=1:length(file_names),.combine=rbind) %dopar%
{
  df=read.delim(file_names[i],header = T)  
  res=colMax(as.data.frame(prop.table(as.matrix(df),margin = 2)*100))
  res=t(as.data.frame(res))
  write.table(res,file=paste(temp_path,file_names[i],sep=""),sep="\t",quote=F,row.names=F,col.names = T)
}}

# Set the frame_counting_matrix directory as current workinf directory and run the function 

setwd("")
frame_percentage(file_names,"path/to/output_dir/")


#Set the output_dir from above step as current directory and run th the followinf command  

setwd("")

listtxt=dir(pattern = "*.txt")

ldf = lapply(listtxt, function(x) {
  dat = read.delim(x, header=T)
  names(dat) = c('27mers', '28mers', '29mers','30mers','31mers','32mers')
  dat$Transcript = x
  return(dat)
})

df=do.call("rbind",ldf)