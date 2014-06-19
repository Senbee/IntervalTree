#############
## Phase 2 ##
#############

Sys.getpid()


### Let's now work on a 'reshaped' annotation object

library(Rcpp)
library(GenomicAlignments)
library(IRanges)


load("TxDb_hg19_UCSC.RData")
load("chr22_list.RData")

length(txdb)
exbytx_list <- split(txdb@unlistData,seqnames(txdb@unlistData))

setClass('IntervalForest_Seq',representation( ptr = "externalptr" ))

setMethod( "initialize", "IntervalForest_Seq", function(.Object, ...) {
    ##.Object@ptr <- makeTree(...)
    .Object@ptr <- makeForest(...)
    .Object
})

Sys.setenv("PKG_CXXFLAGS"="-std=c++11")

## For debugging and optimization 0
## system("cp /usr/local/lib/R/etc/Makeconf.debug cp /usr/local/lib/R/etc/Makeconf")
sourceCpp('IntervalTree.cpp',rebuild=TRUE)

system.time(iForest <- new("IntervalForest_Seq",exbytx_list)) ## ~25 secs


##fl <- BamFile("chr22.bam",asMates=TRUE,qnameSuffixStart="/")
fl1k <- BamFile("chr22.bam",asMates=TRUE,qnameSuffixStart="/",yieldSize=100000)
bam.param <- ScanBamParam(flag=scanBamFlag(isPaired = TRUE, isUnmappedQuery = FALSE, 
                              hasUnmappedMate = FALSE,
                              isNotPassingQualityControls = FALSE),
                          simpleCigar=FALSE,
                          what=scanBamWhat()[c(1,3:5,8,14:15)])
## ## bamWhat(bam.param)

## system.time(chr22_1k <- scanBam(fl1k,param=bam.param)) ## ~ 1 secs
## chr22_1k <- chr22_1k[[1]]
## system.time(chr22 <- scanBam(fl,param=bam.param)) ## ~ 18 secs
## chr22 <- chr22[[1]]


## length(chr22$rname) ## 2328538 reads
## table(chr22$mate_status)
## chr22 <- lapply(chr22,function(x) x[chr22$mate_status == "mated"]) ## can we remove these BEFORE, in scanBam?
## summary(which(duplicated(chr22$groupid)) %% 2) ## all even -> mates are consecutives!

## save(chr22,file="chr22_list.RData")

system.time(res <- getOverlaps(iForest,chr22_1k)) ## ~1.6 secs

