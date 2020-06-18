library(seqinr)
# Loading library seqinr
library(seqinr)

# Loading your input file
fa <- read.fasta("input.fa")

# Getting sequences and gene names from fa object
genes <- sapply(strsplit(names(fa), "\\."), function(v) {return(v[1])})
sequences <- sapply(fa, c2s)

# Extracting longest transcript for each gene
filtered_seq <- tapply(sequences, genes, function(v) {return(v[which(nchar(v)==max(nchar(v)))])})

# Writing output file
write.fasta(filtered_seq , names(filtered_seq), file="output.fasta")