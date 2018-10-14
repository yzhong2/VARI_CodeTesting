source("vcfFunc.R")
library('tidyr')
library('readr')
library('GenomicAlignments')
require("httr")
require("jsonlite")

# read input file
fName <- c('/usr/local/data/data12/yzhong2/test/coding_challenge_final.vcf')
vcfDat <- read_delim(file=fName,delim='\t',col_names = TRUE,comment = "##",progress=FALSE)

#replace the # in the colnames
colnames(vcfDat)[1] <- gsub('#','',colnames(vcfDat)[1])
splName <- colnames(vcfDat)[10]

## process the vcf file, split the multiallelic
dat <- VcfDataProcess(vcfDat=vcfDat,splName=splName)

### specify the out directory and output file name
oFile<-c('coding_challenge_final_Annot.txt')
odir <- c('/usr/local/data/data12/yzhong2/test/')

## annotat the mutation by the Exac using REST service
outFile <- ExAcRest(dat = dat, Odir = odir, fName=oFile)





