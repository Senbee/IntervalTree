Sys.getpid()

library(Rcpp)
library(GenomicAlignments)
library(IRanges)

load(file='Dati.RData')

rm(exbytx) ##let's use only chr4

elementMetadata(exbytx_chr4@unlistData)$exon_name <- paste("ex",seq_len(length(exbytx_chr4@unlistData)),sep="_")
elementMetadata(exbytx_chr4@unlistData)$tx_name <- rep(names(exbytx_chr4),elementLengths(exbytx_chr4))

Sys.setenv("PKG_CXXFLAGS"="-std=c++11")

sourceCpp('IntervalTree.cpp',rebuild=TRUE) 
                                        #,verbose=TRUE,showOutput=TRUE)

setClass('IntervalForest_Seq',representation( ptr = "externalptr" ))
## setMethod( "initialize", "IntervalForest_Seq", function(.Object, ...) {
##     .Object@pointer <- .Call('makeTree', ... )
##     .Object
## })

setMethod( "initialize", "IntervalForest_Seq", function(.Object, ...) {
    .Object@ptr <- makeTree(...)
    .Object
})

## df <- sort(exbytx_chr4@unlistData,ignore.strand=TRUE)
## summary(diff(start(ranges(df))))

iTree <- new("IntervalForest_Seq",exbytx_chr4@unlistData)

id <- 6005L

## read: [26457, 26493] -> OK
## matching interval: [26455, 26667]

## iTree_small <- new("IntervalForest_Seq",exbytx_chr4@unlistData)

tmp <- lapply(s1,function(x) x[6001:6500])


res <- getOverlaps(iTree,tmp)
table(res)

system.time(res <- getOverlaps(iTree,s1)) ## 0.8

