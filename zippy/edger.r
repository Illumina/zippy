###
#Runs edger on our data.
#args
#1: data file
#2: number of samples in group 1
#3: number of samples in group 2
#4: output filename
#5: temp directory
args <- commandArgs(trailingOnly = TRUE)
source("http://bioconductor.org/biocLite.R")
biocLite("edgeR", lib=args[5])
###
library("edgeR", lib.loc=args[5])
#build sample groupings
groups <- c(rep(1,args[2]), rep(2,args[3]))
#read data into edger
print(args[1])
x <- read.delim(args[1], row.names="Symbol", sep=",")
print(groups)
print(colnames(x))
y <- DGEList(counts=x, group=groups)
#filter
keep <- rowSums(cpm(y)>1) >= 1 #very minimal filter step here, requiring a cpm of 1 in at least 1 sample
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y) #performs TMM
design <- model.matrix(~groups)
y <- estimateDisp(y,design)
if(args[2] == 1 || args[3] == 1){
    fit <- glmFit(y,design,dispersion=0.04) #TODO: currently this is somewhat less than the recommended 'human' level of 0.16.  This should be a user-defineable param.
} else {
    fit <- glmFit(y,design)
}
lrt <- glmLRT(fit,coef=2)
write.csv(lrt$table, file=args[4])