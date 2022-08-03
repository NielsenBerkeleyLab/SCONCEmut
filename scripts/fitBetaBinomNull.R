#!/usr/bin/env Rscript
library(extraDistr)
library(stringr)
# Mon 16 May 2022 06:44:59 PM PDT
# script to fit a beta binom distribution to diploid data, by maximizing the loglikelihood of the beta binom
# optimizes loglikelihood by estimating omega (overdispersion)

# P(a_ij | n_ij, f, w) = (n_ij choose a_ij) * Beta(a_ij + fw, n_ij-a_ij + (1-f)w) / Beta(fw, (1-f)w)
# a_ij = # alt reads, for site i, cell j
# n_ij = # total reads, for site i, cell j
# f = allele frequency; but here, using diploid sites not in dbsnp ==> assume is error (ie is not a variable site), so f=1*(1-err) + (1-1)*err = (1-err)
# w = overdispersion parameter (estimated)

# Arg 1: diploid reads matrix (ex dbsnpIntersection/pooledDiploidMajorAllele.llr.dbsnp.notInDbsnp.countSummary.minCells6.2_minAvgReads4.38.readsMat.bed)
# Arg 2: file to get column names from (ex dbsnpIntersection/pooledDiploidMajorAllele.llr.dbsnp.header)
# chr start end dbsnp <cell total num reads> <cell num major allele reads> <cell LLR> [...]
errRate <- 5e-3

args <- commandArgs(trailingOnly=TRUE)
readsTab <- read.table(args[1], header=F)
cols <- system(paste0("head -n 1 ", args[2]), intern=T)
cols <- unlist(str_split(cols, "\t"))
if(cols[length(cols)] == "") {
  cols <- cols[-length(cols)] # remove trailing tab
}
colnames(readsTab) <- cols

# select out numTotalReads matrix
totalReadsTab <- readsTab[,cols[grepl("totalNumReads", cols)]]

# select out numMajorAlleleReads matrix
majorReadsTab <- readsTab[,cols[grepl("numPooledMajorAlleleReads", cols)]]

## subtract to get numMinorAlleleReads mat
#minorReadsTab <- totalReadsTab - majorReadsTab


# allele freq of major allele should be 1 (minus error) ==> care about numMajorReads
f <- 1 - errRate # allele frequency, with error

loglikeFunMax <- function(param) {
  w <- param[1] # omega
  alpha <- f*w
  beta <- (1-f)*w

  # alpha and beta must be nonnegative
  if(alpha < 0 || beta < 0) {
    return(-1e300)
  }
 
  sum(sapply(1:ncol(totalReadsTab), FUN=function(cellIdx) {
    sum(dbbinom(majorReadsTab[,cellIdx], totalReadsTab[,cellIdx], alpha, beta, log=T), na.rm=T)
  }), na.rm=T)
}
par <- 10
res <- optim(par=par, fn=loglikeFunMax, gr=NULL, control=list(fnscale=-1), method="Brent", lower=0, upper=1e6)
sprintf("%.20f", res$par[1])
stop()
print(system.time(res <- optim(par=10, fn=loglikeFunMax, gr=NULL, control=list(fnscale=-1), method="BFGS")))


system.time(optimRes_brent <- lapply(c(1, seq(10, 300, 10)), FUN=function(par) {optim(par=par, fn=loglikeFunMax, gr=NULL, control=list(fnscale=-1), method="Brent", lower=0, upper=1e6)}))
resDf <- do.call(rbind, lapply(optimRes_brent, FUN=function(res) {as.data.frame(t(unlist(res)))}))

# > head(resDf)
#        par     value counts.function counts.gradient convergence
# 1 9.686635 -544547.8              NA              NA           0
# 2 9.686635 -544547.8              NA              NA           0
# 3 9.686635 -544547.8              NA              NA           0
# 4 9.686635 -544547.8              NA              NA           0
# 5 9.686635 -544547.8              NA              NA           0
# 6 9.686635 -544547.8              NA              NA           0
# > summary(resDf)
#       par            value         counts.function counts.gradient  convergence
#  Min.   :9.687   Min.   :-544548   Min.   : NA     Min.   : NA     Min.   :0
#  1st Qu.:9.687   1st Qu.:-544548   1st Qu.: NA     1st Qu.: NA     1st Qu.:0
#  Median :9.687   Median :-544548   Median : NA     Median : NA     Median :0
#  Mean   :9.687   Mean   :-544548   Mean   :NaN     Mean   :NaN     Mean   :0
#  3rd Qu.:9.687   3rd Qu.:-544548   3rd Qu.: NA     3rd Qu.: NA     3rd Qu.:0
#  Max.   :9.687   Max.   :-544548   Max.   : NA     Max.   : NA     Max.   :0
#                                    NA's   :31      NA's   :31
# > sprintf("%.20f", resDf$par[1])
# [1] "9.68663475621146297101"



stop()
##########################

## Arg 1: diploid total reads matrix
## Arg 2: diploid num major allele reads matrix
## Arg 3: output file 

loglikeFunMax <- function(param) {
  w <- param[1] # omega
  alpha <- f*w
  beta <- (1-f)*w
 
  # alpha and beta must be nonnegative
  if(alpha < 0 || beta < 0) {
    return(-1e300)
  }
  loglik <- 0
  for(snpIdx in 1:nrow(totalReadsTab)) {
    for(cellIdx in 1:ncol(totalReadsTab)) {
      total <- totalReadsTab[snpIdx, cellIdx]
      major <- majorReadsTab[snpIdx, cellIdx]
      if(is.na(total) && is.na(major)) {
        next
      }
 
      prob <- dbbinom(major, total, alpha, beta, log=T)
      loglik <- loglik + prob
    }
  }
  return(loglik)
}


print(system.time(res <- optim(10,loglikeFunMax, control=list(fnscale=-1))))
print(system.time(res <- optim(c(10,10),loglikeFunMax, control=list(fnscale=-1))))
#print(res)
sprintf("%.20f", res$par)
alphaBetaTab <- data.frame(c("alpha", "beta"), res$par)
write.table(format(alphaBetaTab, digits=10), file=outputFile, quote=F, sep="=", col.names=F, row.names=F)

# calc ll for all snps and save to file
# TODO

loglikeFunMaxApply <- function(param) {
  alpha <- param[1]
  beta <- param[2]
  # alpha and beta must be nonnegative
  if(alpha < 0 | beta < 0) {
    return(-1e300)
  }
  #loglik <- sum(sapply(1:1000000, FUN=function(snpIdx) {
  loglik <- sum(sapply(1:nrow(totalReadsTab), FUN=function(snpIdx) {
    sum(sapply(1:ncol(totalReadsTab), FUN=function(cellIdx) {
      total <- totalReadsTab[snpIdx, cellIdx]
      major <- majorReadsTab[snpIdx, cellIdx]
      if(is.na(total)) {
        return(0)
      }
      dbbinom(major, total, alpha, beta, log=T)
    }))
  }))
  loglik
}

