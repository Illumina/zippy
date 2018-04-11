###
#Runs deseq2 on our data.
#TODOS:
##We can only normalize against a single sample
##We have chosen a fittype that may not make sense for real data
#args
#1: data directory
#2: control sample
#3: output filename
#source("http://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
###
library("DESeq2")
args <- commandArgs(trailingOnly = TRUE)
#get output from htseq
sampleFiles <- grep("hts.txt",list.files(args[1]),value=TRUE)
#condition is currently based on file name.  We may eventually want to take in a list for the condition
sampleCondition <- sub("(.*).hts.txt","\\1",sampleFiles)
sampleTable <- data.frame(sampleName = sampleFiles, fileName = sampleFiles, condition = sampleCondition)
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = args[1], design= ~ condition)

dds <- dds[ rowSums(counts(dds)) > 1, ] # filter to genes with at least 2 reads
dds <- DESeq(dds, fitType="mean")
dds$condition <- relevel(dds$condition, ref=args[2]) # define the control sample.
res <- results(dds)
write.csv(as.data.frame(res), file=paste(args[3], "deseq.csv", sep="/"))
