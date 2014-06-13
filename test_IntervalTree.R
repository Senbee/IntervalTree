## TODO:
## R code:
##  1) check seqnames consistency between reads & txdb
##
## C++ Code:
##  1) define IntervalForest class to store as many trees as chromosomes
##  2) create search method for IntervalForest (in the class!)


Sys.getpid()

library(Rcpp)
library(GenomicAlignments)
library(IRanges)

load(file='Dati.RData')


exbytx_list <- split(exbytx@unlistData,seqnames(exbytx@unlistData))

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

## id <- 6005L

## read: [26457, 26493] -> OK
## matching interval: [26455, 26667]

## iTree_small <- new("IntervalForest_Seq",exbytx_chr4@unlistData)

tmp <- lapply(s1,function(x) x[6001:6500])

res <- getOverlaps(iTree,tmp)
table(res)

system.time(res <- getOverlaps(iTree,s1)) ## ~0.8







### Let's now work on a 'reshaped' annotation object


source('~/Work/Sweden/RNA-seq/Package/Sequgio/R/reshapeTxDb.R')

reshape(exbytx_chr4)

