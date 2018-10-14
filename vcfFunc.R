## this function is for process vcf data set

VcfDataProcess <- function(vcfDat,splName) {
   count <- 0 
   proDat <- apply(vcfDat, 1, function(r) {
             ## process allete genotype
             count <<- count + 1
             ## this function is to split multiallelic case into separated rows
             alleleMat <- alleGenoType(oneLine=r,cName="ALT")
             repTimes <- nrow(alleleMat)
             ## get the reference genotype
             refMat <-CoordChromPro(oneLine=r,repTimes=repTimes,cName="REF")
             ## get the chromesome ID for allele genotype
             chrMat <-CoordChromPro(oneLine=r,repTimes=repTimes,cName="CHROM") 
             ## get the starting position for allele genotype
             posMat<- CoordChromPro(oneLine=r,repTimes=repTimes,cName="POS")
             ## get the Cigar value
             cigarMat <- colInfoTag(oneLine=r,cName='INFO',tagName='CIGAR')
             ## get the Allele frequency
             AFMat <- colInfoTag(oneLine=r,cName='INFO',tagName='AF')
             ## get the read number for ref and alleles
             RDpMat<-readAlleleDepth(oneLine=r,cName='FORMAT', splName=splName) 
             resultMat <- cbind(chrMat,posMat,refMat,alleleMat,cigarMat,
                                AFMat,RDpMat)
             ## this step is for cigar value splitting for multiple muations
             ## happen in one cigar. For example cigar value : 2X5D2X, we need to split this row into 5 rows
             ## with each row for one mutation.
             cigarMat <- apply(resultMat,1,CigarValuePro)
             cigarMat <- cigarMat[[1]][[1]]
             colnames(cigarMat) <- colnames(resultMat)
             return (cigarMat)
             })
  mutdat <- do.call('rbind',proDat)
  ## this step I will add the mutation type
  mutdat <- VariantType(mutDat=mutdat, refName='REF', allName='ALT')
  rownames(mutdat) <- NULL
  return (mutdat)
}

### this is for vriant type annotation
VariantType <- function(mutDat, refName, allName) {
  ## this is for the insertion
  rlen <- nchar(mutDat[,refName])
  qlen <- nchar(mutDat[,allName])
  mutType <- vector(mode = 'character',length=nrow(mutDat))
  index <- rlen == qlen
  mutType[index] <- 'SNP'
  index <- rlen > qlen
  mutType[index] <- 'DEL'
  index <- rlen < qlen
  mutType[index] <- 'INS'
  mutDat <- cbind(mutDat,cbind(mutType))
  return (mutDat)
}





### this function is for the cigar value process
### need to split one row into multiple rows if cigar value contain multiple
### mutation

CigarValuePro <- function(oneLine,cName='CIGAR') {
   cigar <- oneLine[cName]
   cigarOps <- explodeCigarOps(cigar, ops=CIGAR_OPS)[[1]]
   cigarLen <- explodeCigarOpLengths(cigar, ops=CIGAR_OPS)[[1]]
   cigarTbl <- cigarOpTable(cigar)
   if(sum(cigarTbl[,c('N','S','H','P', '=')]) != 0) {
     stop('need to process this case for cigar value')
   }
   result <- PostionCigar(oneLine=oneLine,refSeq=oneLine['REF'],
                          qurySeq=oneLine['ALT'], cigarOpe=cigarOps,
                          cigarLen=cigarLen)
   return (result)
}

## this function I will calculate the mutation position by the cigar value 
PostionCigar <- function(oneLine,refSeq, qurySeq,cigarOpe,cigarLen) {
   oneLineMat <- rbind(oneLine)
   colnames(oneLineMat) <- names(oneLine)
   refbPos <- 0 
   refePos <- 0 
   refbVect <- c()
   refeVect <- c()
   qurybPos <- 0 
   quryePos <- 0 
   qurybVec <- c()
   quryeVec <- c()
   refGeno <- c()
   quryGeno <- c()
   for(i in 1 : length(cigarOpe)) {
      opType <- cigarOpe[i]
      opNum  <- cigarLen[i]
      if(opType == 'X') {
         ### process the mismatch
         dispath <- rep(1,times=opNum)
         for(m in 1 : length(dispath)) {
            refbPos <- refbPos + dispath[m] 
            refePos <- refbPos 
            refbVect <- c(refbVect,refbPos)
            refeVect <- c(refeVect,refePos) 
            qurybPos <- qurybPos + dispath[m] 
            quryePos <- qurybPos 
            qurybVec <- c(qurybVec,qurybPos)
            quryeVec <- c(quryeVec,quryePos)
            refGeno <- c(refGeno,substr(refSeq, refbPos, refePos))
            quryGeno <- c(quryGeno,substr(qurySeq,qurybPos,quryePos))
          }
      }
      else if(opType == 'M') {
         refbPos <- refbPos + opNum 
         refePos <- refbPos 
         qurybPos <- qurybPos + opNum 
         quryePos <- qurybPos 
      }
      ## this is for deletion , so reference should have ref geno, query 
      ## geno should be short, like "chr2	61414639	GCAAT	G 1M4D deletion"
      else if(opType == 'D'){ 
         refePos <- refePos + opNum    
         refGeno <- c(refGeno,substr(refSeq, refbPos, refePos)) 
         quryGeno <- c(quryGeno,substr(qurySeq,qurybPos,qurybPos))
         qlen <- length(quryGeno)
         relen <- length(refGeno)
         if(substr(quryGeno[qlen],1,1) != substr(refGeno[relen],1,1)) {
            stop('deletion genotype is wrong')
         }
         refbVect <- c(refbVect,refbPos)
         refeVect <- c(refeVect,refePos)
         refbPos <- refePos + 1
         refePos <- refbPos
      }
      ## this step is for the insertion case, so reference genotype should 
      ## be shorter but qury should contain insertion sequence like : chr2	209223123	G	GTTT	1M3I insertion
      else if(opType == 'I'){
         quryePos <- quryePos + opNum 
         quryGeno <- c(quryGeno,substr(qurySeq, qurybPos,quryePos))
         quryePos <- quryePos + 1
         qurybPos <- quryePos
         ## this is for the reference geno
         refGeno <- c(refGeno, substr(refSeq,refbPos,refePos))
         qlen <- length(quryGeno)
         relen <- length(refGeno)
         if(substr(quryGeno[qlen],1,1) != substr(refGeno[relen],1,1)) {
            stop('insertion genotype is wrong')
         }
         refbVect <- c(refbVect,refbPos)
         refeVect <- c(refeVect,refePos)
      }
      else {
         stop('we dont have this genotype')
      }
   }
   startpos <- as.integer(oneLineMat[,'POS'])
   if(length(refbVect) > 1) {
      rowIndex <- rep(1, times=length(refbVect))
      oneLineMat <- oneLineMat[rowIndex,]
   }
   oneLineMat[,'POS'] <- as.character(refbVect + startpos - 1)
   oneLineMat[,'REF'] <- refGeno
   oneLineMat[,'ALT'] <- quryGeno
   return (list(oneLineMat))
}


### this function is to split the allele genotype for the multi-allele cases
alleGenoType <- function(oneLine,cName) {
   allSeqList <- strsplit(oneLine[cName],',')[[1]]
   alleMat <- cbind(allSeqList)
   colnames(alleMat) <- cName
   return (alleMat)
}

### this function is for the chrom and coordinate process
CoordChromPro <- function(oneLine, repTimes, cName) {
   eleM <- cbind(rep(oneLine[cName], times=repTimes))
   colnames(eleM) <- cName
   return (eleM)
}

## this function is usd for process INFO column in vcf file
colInfoTag <- function(oneLine,cName,tagName) {
  info <- oneLine[cName]
  index <- grep(tagName, info)
  if(length(index) == 0) {
     stop(tagName,' not exist')
  }
  pattern <- paste(tagName,'=',sep='')
  tags <- gsub(pattern,'',grep(tagName,strsplit(info,';')[[1]],value=TRUE))
  ## this split the tags further for multi-allele case
  tagsMat <- cbind(strsplit(tags,',')[[1]])
  colnames(tagsMat) <- tagName 
  return (tagsMat)
}

### this function is for extracting the read depth and reads support allele
### and reads supporting the reference
readAlleleDepth <- function(oneLine, cName, splName) {
  ForMat <- oneLine[cName]
  splInfo <- oneLine[splName]
  ## get the tag names   
  names <- strsplit(ForMat,':')[[1]]
  readNum <- strsplit(splInfo,':')[[1]]
  names(readNum) <- names
  DpRead <- readNum['DP']
  RoRead <- readNum['RO']
  AoRead <- readNum['AO'] 
  ### this step is to split the AoRead for multi-allele case
  AoRead <- cbind(strsplit(AoRead,',')[[1]])
  colnames(AoRead) <- 'alleleReadDepth'
  mode(AoRead) <- 'numeric'
  repTimes <- nrow(AoRead)
  DpRead <- cbind(rep(DpRead, times=repTimes))
  colnames(DpRead) <-'TotalReadDepth'
  mode(DpRead) <- 'numeric'
  RoRead <- cbind(rep(RoRead, times=repTimes))
  colnames(RoRead) <- 'RefReadDepth' 
  mode(RoRead) <- 'numeric'
  alleleProportion <- AoRead/DpRead
  colnames(alleleProportion) <- c('alleleProportion')
  refProportion <- RoRead/DpRead
  colnames(refProportion) <- c('refProportion')
  resultMat <- cbind(DpRead,RoRead,refProportion,AoRead,alleleProportion)
  return (resultMat)
}

### this function is about the Exac database access
### http://exac.broadinstitute.org/variant/X-153296070-A-AG
ExAcRest <- function(dat, url='http://exac.hms.harvard.edu/rest/variant/variant/', Odir, fName) {
  cName <- c('SYMBOL','major_consequence','SOMATIC','HGVSc','HGVSp','LoF_info',
             "LoF_filter","LoF_flags",'Existing_variation','PolyPhen','SIFT')
  annotDat  <- apply(dat, 1, function(r) {
            chrId <- gsub('chr','',r['CHROM'])
            pos <- r['POS']
            ref <- r['REF']
            alt <- r['ALT']
            chrs <- paste(url,paste(chrId,pos,ref,alt,sep='-'),sep='')
            oneResult <- content(GET(chrs),"text",encoding ="UTF-8")
            json <- fromJSON(oneResult, flatten = TRUE)
            jsonDat <- json$vep_annotations
            #diffcol <- setdiff(cName, colnames(jsonDat))
            if(length(jsonDat) == 0) {
              nullDat<-as.data.frame(rbind(rep('NoExac',times=length(cName))))
              colnames(nullDat) <- cName
              return (nullDat)
            }
            else {
              seleJson <- jsonDat[,cName]
              seleJson <- selectionDeleter(dat=seleJson)
            }
            return (seleJson)
            })
   annotDat <- do.call('rbind', annotDat)
   rownames(annotDat) <- NULL 
   wholeDat <- cbind(dat, annotDat)
   ## this I will write dat file
   fName <- file.path(Odir, fName)

   write.table(wholeDat,file =fName, quote = FALSE, sep = "\t",
              row.names = FALSE, col.names = TRUE )
   return (fName)
}

### this function is for selecting the multiple mutation function impacts
### selection will be based on polyphen score
selectionDeleter <- function(dat) {
   if(nrow(dat) == 1) {
     return (dat)
   }
   else {
     index <- which(dat[,'PolyPhen'] !="")
     if(length(index) == 0) {
       return (dat[1,,drop=FALSE])
     }
     else {
       tempDat <- dat[index,,drop=FALSE]
       score <- apply(tempDat,1, function(r){
       matches <- parse_number(r['PolyPhen']) 
       return (matches)
       })
       mIndex <- which(score==max(score))[1]
       return (tempDat[mIndex,,drop=FALSE])
     }
   }
}






