
    data = read.table("GC_content.dat", header=F, row.names=1)
    pdf("GC_content.dat.hist.pdf")
    hist(data[,1], br=100)
    message("

mean: ", sprintf("%.2f", mean(data[,1])), ", median: ", median(data[,1]), "

")
    dev.off()

