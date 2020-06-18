################################## fOverlaps ###############################################
temp_df=data.frame(start=11150,end=11170,strand="-")

last_overlapping_read=fucntion()
{

  for(i in 1:dim(df)[1])
    {
      x=data.table(seqnames=temp$seqnames,start=temp$start,end=temp$end,strand=temp$strand)
      y=data.table(seqnames=df$seqnames,start=df$start,end=df$end,strand=df$strand)
      setkey(x,seqnames,start,end)
      setkey(y,seqnames,start,end)
      temp=foverlaps(x, y, type="any",mult="last")
      res=temp
    }
  return() 
}
print(res)

################################### findOverlaps ##############################################

start_position=1151
end_position=1170

temp_df=data.frame(start=11150,end=11170)
for(i in 1:dim(df)[1])
{
  temp_df=df[subjectHits(findOverlaps(IRanges(temp_df$start[dim(temp_df)[1]], temp_df$end[dim(temp_df)[1]]), IRanges(df$start, df$end))),]
  res=temp_df
}
print(res)














