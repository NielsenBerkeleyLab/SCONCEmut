library(ape)
library(plyr)
library(reshape2)
library(stringr)
library(ggplot2)
library(cowplot)
library(scales)

# Thu 10 Feb 2022 03:02:54 PM PST
# script of shared tree reading functions and constants (newick strings, branch lengths)

###########################################
# Tree manipulations for trees read from newick strings
###########################################
# for sp21, Node 127 and Node 128 correspond to labels 1 and 2 in newick file. Node 127 is cellNum 127, Node 128 is cellNum 126; all the way to Node 254 == celNum 0 == label 128
# node 127 == label 1 == cellnum 127
# node 128 == label 2 == cellnum 126
# node 254 == label 128 == cellnum 0
# assumes sp22 is the same. or if not, doesn't really matter since all equal height
convertCellNumToNewickLabel <- function(cellFilename) { # ex simu_cancer_cell_0.hg19_lite.depth
  cellNum <- as.numeric(str_extract(str_extract(cellFilename, "cell_[0-9]+"), "[0-9]+"))
  label <- 128 - cellNum
  #label <- 5 - cellNum
  label
}

# given mrcaMat, distMat, nodeDepths, and branchLength, calculates all possible tree branch combinations
# left | right | t1 | t2 | t3
calcTreeBranches <- function(mrcaMat, distMat, nodeDepths, branchLength) {
  do.call(rbind, lapply(1:(ncol(mrcaMat)), FUN=function(left) {
    do.call(rbind, lapply(1:ncol(mrcaMat), FUN=function(right) {
      currMrca <- mrcaMat[left, right]
      #t1 <- nodeDepths[currMrca] + branchLength # add for time before first split
      t1 <- nodeDepths[currMrca]# + branchLength # add for time before first split Thu 30 Jun 2022 11:30:16 AM PDT I think this might be wrong?
      #t1 <- nodeDepths[currMrca] + 1# + branchLength # add for time before first split Thu 30 Jun 2022 11:30:16 AM PDT I think this might be wrong?
      t2 <- distMat[left, currMrca]
      t3 <- distMat[right, currMrca]
      results <- data.frame(left=left, right=right, t1=t1, t2=t2, t3=t3, t2_t3=t2+t3)
      results
    }))
  }))
}

###########################################
# File reading functions
###########################################
# reads .hmm files to get tree branch estimates and generates true tree branches; returns a matrix with
# left (newick label) | right | variable (branch number) | cell0 (filename) | cell1 | inferred | true
getTreeBranches <- function(newickString, branchLength, hmmFile, forceRecalc=F) {
  outputFile <- paste0(hmmFile, ".treeBranches.txt")
  if(!forceRecalc && file.exists(outputFile)) {
    merged <- read.table(outputFile, sep="\t", header=T, stringsAsFactors=F)
    return(merged)
  }
  treeObj <- read.tree(text=newickString)
  mrcaMat <- mrca(treeObj)
  distMat <- dist.nodes(treeObj)
  nodeDepths <- node.depth.edgelength(treeObj)

  # left | right | t1 | t2 | t3
  # left/right correspond to node labels in the newick file
  treeBranches <- calcTreeBranches(mrcaMat, distMat, nodeDepths, branchLength)
  treeBranches_m <- melt(treeBranches, id.vars=c("left", "right"))

  # cell0 | cell1 | t1
  # cell0 | cell1 | t2
  # cell0 | cell1 | t3
  # ...
  # need to separate cell name reading from param reading in case doing nearest n
  hmmLines <- system(paste0("/bin/bash -c ", shQuote(sprintf("grep \"^HMM\" %s | grep \",\" | sed -e 's/HMM\\s\\+//' -e 's/(//' -e 's/)//' -e 's/,/\t/' -e 's/://' -e 's/\\s\\+/\t/'", hmmFile))), intern=T)
  if(length(hmmLines) == 0) {
    return(NULL)
  }
  analyzedHmmNames <- read.table(text=hmmLines)
  colnames(analyzedHmmNames) <- c("hmmIdx", "cell0", "cell1")
  paramIdxToExtract <- do.call(c, lapply(analyzedHmmNames$hmmIdx, FUN=function(i) {c(3*i, (3*i)+1, (3*i)+2)})) # don't extract hmms we skipped
  paramIdxToExtract <- paramIdxToExtract + 1 # paramsToEst is 0 indexed, R vector is 1 indexed
  unlabeledHmmParamsRaw <- read.table(text=system(paste0("/bin/bash -c ", shQuote(sprintf("sed '1,/FINAL HMM/d' %s | tac | sed '/paramsToEst/q' | tac | sed '/fixedParams/q' | tail -n +2 | head -n -2", hmmFile))), intern=T))
  hmmParamsRaw <- data.frame(cell0=do.call(c, lapply(analyzedHmmNames$cell0, rep, 3)), cell1=do.call(c, lapply(analyzedHmmNames$cell1, rep, 3)), value=unlabeledHmmParamsRaw[paramIdxToExtract,])
  
  hmmParamsRaw$variable <- c("t1", "t2", "t3")

  # cell0 | cell1 | t1 | t2 | t3
  # cell0/cell1 correspond to indv cell filenames
  hmmParams_wide <- dcast(hmmParamsRaw, cell0 + cell1 ~ variable)
  hmmParams_wide$t2_t3 <- hmmParams_wide$t2 + hmmParams_wide$t3
  #hmmParams_wide$total <- hmmParams_wide$t1 + hmmParams_wide$t2 + hmmParams_wide$t3
  ##hmmParams_wide$total <- hmmParams_wide$t1 + max(hmmParams_wide$t2,hmmParams_wide$t3)
  ##hmmParams_wide$total <- hmmParams_wide$t2 + hmmParams_wide$t3
  #hmmParams_wide$t1_sc <- (hmmParams_wide$t1) / (hmmParams_wide$total)
  #hmmParams_wide$t2_sc <- (hmmParams_wide$t2) / (hmmParams_wide$total)
  #hmmParams_wide$t3_sc <- (hmmParams_wide$t3) / (hmmParams_wide$total)
  #hmmParams_wide$t2_t3_sc <- (hmmParams_wide$t2 + hmmParams_wide$t3) / (hmmParams_wide$total)
  hmmParams_wide_m <- melt(hmmParams_wide, id.vars=c("cell0", "cell1"))

  hmmParams_wide_m$left <- sapply(hmmParams_wide_m$cell0, convertCellNumToNewickLabel)
  hmmParams_wide_m$right <- sapply(hmmParams_wide_m$cell1, convertCellNumToNewickLabel)

  merged <- merge(hmmParams_wide_m, treeBranches_m, by=c("left", "right", "variable"))
  colnames(merged) <- c("left", "right", "variable", "cell0", "cell1", "inferred", "true")

  write.table(merged, file=outputFile, quote=F, sep="\t", col.names=T, row.names=F)

  merged
}

# based on non zero alleles directly from simulations
getTruePerBranchMutCounts <- function(newickString, branchLength, paramSet, forceRecalc=F, uniqOnly=T, numCellsList=NULL) {
  uniqFilt <- ifelse(uniqOnly, "", ".all")
  numCellsFilt <- ifelse(is.null(numCellsList), "", paste0("_c", paste0(numCellsList, collapse="-c")))
  outputFile <- paste0(dataDir, paramSet, "/truePerBranchMutCounts", uniqFilt, numCellsFilt)
  if(!forceRecalc && file.exists(paste0(outputFile, "_indv.txt")) && file.exists(paste0(outputFile, "_pairs.txt"))) {
    indvDat <- read.table(paste0(outputFile, "_indv.txt"), sep="\t", header=T)
    pairDat <- read.table(paste0(outputFile, "_pairs.txt"), sep="\t", header=T)
    return(list(indv=indvDat, paired=pairDat))
  }
  # left | right | t1 | t2 | t3
  # left/right correspond to node labels in the newick file
  treeObj <- read.tree(text=newickString)
  mrcaMat <- mrca(treeObj)
  distMat <- dist.nodes(treeObj)
  nodeDepths <- node.depth.edgelength(treeObj)
  treeBranches <- calcTreeBranches(mrcaMat, distMat, nodeDepths, branchLength)
  treeBranches_m <- melt(treeBranches, id.vars=c("left", "right"))

  trueSnpAllelesFileList <- system(paste0("find ", dataDir, paramSet, " -maxdepth 1 -name \"true_snpAlleles_cancer_cell_*\" | sort -V"), intern=T)
  alleleCountsList <- lapply(trueSnpAllelesFileList, FUN=function(f) {
    cellName <- basename(f)
    tab <- read.table(f, sep="\t", header=F)
    colnames(tab) <- c("snpIdx", paste0(cellName, "_ancAlleles"), paste0(cellName, "_derAlleles"))
    tab[,2:3]
  })
  cellNames <- sapply(trueSnpAllelesFileList, basename)
  names(alleleCountsList) <- cellNames

  # for each individual cell
  allIndMutCounts <- do.call(rbind, lapply(1:length(cellNames), FUN=function(cellIdx) {
    cell <- cellNames[cellIdx]
    numMuts <- sum(alleleCountsList[[cellIdx]][,2] != 0)
    left <- convertCellNumToNewickLabel(cell)
    #branchLength <- nodeDepths[left] + branchLength
    branchLength <- nodeDepths[left]# + branchLength Thu 30 Jun 2022 03:01:21 PM PDT I think this may be wrong?
    data.frame(cell=cell, t_muts=numMuts, left=left, true=branchLength)
  }))

  mutCounts <- NULL
  if(is.null(numCellsList)) {
    # for each pair
    if(uniqOnly) {
      cell0IdxRange <- 1:(length(cellNames) - 1)
    } else  {
      cell0IdxRange <- 1:length(cellNames)
    }
    mutCounts <- do.call(rbind, lapply(cell0IdxRange, FUN=function(cell0Idx) {
      if(uniqOnly) {
        cell1IdxRange <- (cell0Idx + 1):length(cellNames)
      } else  {
        cell1IdxRange <- 1:length(cellNames)
      }
      do.call(rbind, lapply(cell1IdxRange, FUN=function(cell1Idx) {
        # calc number germline muts (t0), shared muts (t1), and number muts in only cell (t2/t3)
        # germline: [*, 0],  [*,0] (note, these are later ignored)
        numGermline <- sum(alleleCountsList[[cell0Idx]][, 2] == 0 & alleleCountsList[[cell1Idx]][, 2] == 0)
        # t1:       [*, !0], [*, !0]
        numT1 <- sum(alleleCountsList[[cell0Idx]][, 2] != 0 & alleleCountsList[[cell1Idx]][, 2] != 0)
        # t2:       [*, !0], [*, 0]
        numT2 <- sum(alleleCountsList[[cell0Idx]][, 2] != 0 & alleleCountsList[[cell1Idx]][, 2] == 0)
        # t3:       [*, 0],  [*, !0]
        numT3 <- sum(alleleCountsList[[cell0Idx]][, 2] == 0 & alleleCountsList[[cell1Idx]][, 2] != 0)
        data.frame(cell0=cellNames[cell0Idx], cell1=cellNames[cell1Idx], germline=numGermline, t1=numT1, t2=numT2, t3=numT3)
      }))
    }))
  } else {
    mutCounts <- do.call(rbind, lapply(numCellsList, FUN=function(numCells) {
      currDepths <- read.table(paste0(dataDir, paramSet, "/tumor_depths_", numCells))
      cellNames <- basename(currDepths$V1)
      cellNames <- gsub(".hg19_lite.depth", "", gsub("simu_readDepth", "true_snpAlleles", cellNames))
      cellNums <- as.numeric(str_extract(str_extract(cellNames, "cell_[0-9]+"), "[0-9]+"))
      alleleCountsList <- lapply(cellNames, FUN=function(cellName) {
        tab <- read.table(paste0(dataDir, paramSet, "/", cellName), sep="\t", header=F)
        colnames(tab) <- c("snpIdx", paste0(cellName, "_ancAlleles"), paste0(cellName, "_derAlleles"))
        tab[,2:3]
      })
      do.call(rbind, lapply(1:(length(cellNames)-1), FUN=function(cell0Idx) {
        do.call(rbind, lapply((cell0Idx+1):length(cellNames), FUN=function(cell1Idx) {
          # calc number germline muts (t0), shared muts (t1), and number muts in only cell (t2/t3)
          # germline: [*, 0],  [*,0] (note, these are later ignored)
          numGermline <- sum(alleleCountsList[[cell0Idx]][, 2] == 0 & alleleCountsList[[cell1Idx]][, 2] == 0)
          # t1:       [*, !0], [*, !0]
          numT1 <- sum(alleleCountsList[[cell0Idx]][, 2] != 0 & alleleCountsList[[cell1Idx]][, 2] != 0)
          # t2:       [*, !0], [*, 0]
          numT2 <- sum(alleleCountsList[[cell0Idx]][, 2] != 0 & alleleCountsList[[cell1Idx]][, 2] == 0)
          # t3:       [*, 0],  [*, !0]
          numT3 <- sum(alleleCountsList[[cell0Idx]][, 2] == 0 & alleleCountsList[[cell1Idx]][, 2] != 0)
          data.frame(cell0=cellNames[cell0Idx], cell1=cellNames[cell1Idx], germline=numGermline, t1=numT1, t2=numT2, t3=numT3)
        }))
      }))
    }))
  }
  rownames(mutCounts) <- NULL

  mutCounts_m <- melt(mutCounts, id.vars=c("cell0", "cell1"))
  mutCounts_m$left <-  sapply(mutCounts_m$cell0, convertCellNumToNewickLabel)
  mutCounts_m$right <- sapply(mutCounts_m$cell1, convertCellNumToNewickLabel)

  merged <- merge(mutCounts_m, treeBranches_m, by=c("left", "right", "variable"))
  colnames(merged) <- c("left", "right", "variable", "cell0", "cell1", "inferred", "true")

  write.table(allIndMutCounts, file=paste0(outputFile, "_indv.txt"), quote=F, sep="\t", col.names=T, row.names=F)
  write.table(merged, file=paste0(outputFile, "_pairs.txt"), quote=F, sep="\t", col.names=T, row.names=F)

  list(indv=allIndMutCounts, paired=merged)
}

# based on non zero read counts directly from simulations
getObservedPerBranchMutCounts <- function(newickString, branchLength, paramSet, mutFilt, forceRecalc=F, uniqOnly=T, numCellsList=numCellsList) {
  mutFilt <- gsub("^_", ".", mutFilt)
  uniqFilt <- ifelse(uniqOnly, "", ".all")
  numCellsFilt <- ifelse(is.null(numCellsList), "", paste0("_c", paste0(numCellsList, collapse="-c")))
  outputFile <- paste0(dataDir, paramSet, "/obsPerBranchMutCounts", mutFilt, uniqFilt, numCellsFilt)
  if(!forceRecalc && file.exists(paste0(outputFile, "_indv.txt")) && file.exists(paste0(outputFile, "_pairs.txt"))) {
    indvDat <- read.table(paste0(outputFile, "_indv.txt"), sep="\t", header=T)
    pairDat <- read.table(paste0(outputFile, "_pairs.txt"), sep="\t", header=T)
    return(list(indv=indvDat, paired=pairDat))
  }
  # left | right | t1 | t2 | t3
  # left/right correspond to node labels in the newick file
  treeObj <- read.tree(text=newickString)
  mrcaMat <- mrca(treeObj)
  distMat <- dist.nodes(treeObj)
  nodeDepths <- node.depth.edgelength(treeObj)
  treeBranches <- calcTreeBranches(mrcaMat, distMat, nodeDepths, branchLength)
  treeBranches_m <- melt(treeBranches, id.vars=c("left", "right"))

  obsSnpReadsFileList <- system(paste0("find ", dataDir, paramSet, " -maxdepth 1 -regex \".*/simu_snpAlleles_cancer_cell_[0-9]+", mutFilt, ".hg19_lite.bed\" | sort -V"), intern=T)
  readsCountsList <- lapply(obsSnpReadsFileList, FUN=function(f) {
    cellName <- basename(f)
    tab <- read.table(f, sep="\t", header=F)
    colnames(tab) <- c("chr", "start", "end", paste0(cellName, "_ancReads"), paste0(cellName, "_derReads"))
    tab$key <- paste0(tab$chr, ":", tab$start, "-", tab$end)
    tab
  })
  cellNames <- sapply(obsSnpReadsFileList, basename)
  names(readsCountsList) <- cellNames

  # for each individual cell
  allIndMutCounts <- do.call(rbind, lapply(1:length(cellNames), FUN=function(cellIdx) {
    cell <- cellNames[cellIdx]
    numMuts <- sum(readsCountsList[[cellIdx]][,2] != 0)
    left <- convertCellNumToNewickLabel(cell)
    #branchLength <- nodeDepths[left] + branchLength
    branchLength <- nodeDepths[left]# + branchLength Thu 30 Jun 2022 03:01:21 PM PDT I think this may be wrong?
    data.frame(cell=cell, t_muts=numMuts, left=left, true=branchLength)
  }))

  mutCounts <- NULL
  if(is.null(numCellsList)) {
    # for each pair
    if(uniqOnly) {
      cell0IdxBound <- length(cellNames) - 1
    } else  {
      cell0IdxBound <- length(cellNames)
    }
    mutCounts <- do.call(rbind, lapply(1:cell0IdxBound, FUN=function(cell0Idx) {
      if(uniqOnly) {
        cell1IdxBound <- cell0Idx + 1
      } else  {
        cell1IdxBound <- 1
      }
      do.call(rbind, lapply(cell1IdxBound:length(cellNames), FUN=function(cell1Idx) {
        # calc number germline muts (t0), shared muts (t1), and number muts in only cell (t2/t3)
        derReadsCol <- which(grepl("derReads", colnames(readsCountsList[[cell0Idx]]))) # assumes same for both cells
        # get muts with data in both cells
        inBoth <- list()
        inBoth[[1]] <- readsCountsList[[cell0Idx]][readsCountsList[[cell0Idx]]$key %in% readsCountsList[[cell1Idx]]$key,]
        inBoth[[2]] <- readsCountsList[[cell1Idx]][readsCountsList[[cell1Idx]]$key %in% readsCountsList[[cell0Idx]]$key,]

        # get muts with data in only one cell
        onlyCell0 <- readsCountsList[[cell0Idx]][!(readsCountsList[[cell0Idx]]$key %in% readsCountsList[[cell1Idx]]$key),]
        onlyCell1 <- readsCountsList[[cell1Idx]][!(readsCountsList[[cell1Idx]]$key %in% readsCountsList[[cell0Idx]]$key),]

        # to be in t1, must have data in both cells and have !0 in der alleles
        numT1 <- sum(inBoth[[1]][, derReadsCol] != 0 & inBoth[[2]][, derReadsCol] != 0)
        # to be in t2, must have (data in only cell0 and !0 in der alleles) || (data in both cells and !0 in cell0 der but 0 in cell1 der)
        numT2 <- sum(onlyCell0[, derReadsCol] != 0) + sum(inBoth[[1]][, derReadsCol] != 0 & inBoth[[2]][, derReadsCol] == 0)
        # to be in t3, must have (data in only cell3 and !0 in der alleles) || (data in both cells and !0 in cell1 der but 0 in cell0 der)
        numT3 <- sum(onlyCell1[, derReadsCol] != 0) + sum(inBoth[[1]][, derReadsCol] == 0 & inBoth[[2]][, derReadsCol] != 0)

        data.frame(cell0=cellNames[cell0Idx], cell1=cellNames[cell1Idx], t1=numT1, t2=numT2, t3=numT3)
      }))
    }))
  } else {
    mutCounts <- do.call(rbind, lapply(numCellsList, FUN=function(numCells) {
      currReadsCounts <- read.table(paste0(dataDir, paramSet, "/tumor_depths_", numCells))
      cellNames <- basename(currReadsCounts$V1)
      cellNames <- gsub(".hg19_lite.depth", paste0(mutFilt, ".hg19_lite.bed"), gsub("readDepth", "snpAlleles", cellNames))
      readsCountsList <- lapply(cellNames, FUN=function(cellName) {
        tab <- read.table(paste0(dataDir, paramSet, "/", cellName), sep="\t", header=F)
        colnames(tab) <- c("chr", "start", "end", paste0(cellName, "_ancReads"), paste0(cellName, "_derReads"))
        tab$key <- paste0(tab$chr, ":", tab$start, "-", tab$end)
        tab
      })
      do.call(rbind, lapply(1:(length(cellNames)-1), FUN=function(cell0Idx) {
        do.call(rbind, lapply((cell0Idx+1):length(cellNames), FUN=function(cell1Idx) {
          # calc number germline muts (t0), shared muts (t1), and number muts in only cell (t2/t3)
          derReadsCol <- which(grepl("derReads", colnames(readsCountsList[[cell0Idx]]))) # assumes same for both cells
          # get muts with data in both cells
          inBoth <- list()
          inBoth[[1]] <- readsCountsList[[cell0Idx]][readsCountsList[[cell0Idx]]$key %in% readsCountsList[[cell1Idx]]$key,]
          inBoth[[2]] <- readsCountsList[[cell1Idx]][readsCountsList[[cell1Idx]]$key %in% readsCountsList[[cell0Idx]]$key,]

          # get muts with data in only one cell
          onlyCell0 <- readsCountsList[[cell0Idx]][!(readsCountsList[[cell0Idx]]$key %in% readsCountsList[[cell1Idx]]$key),]
          onlyCell1 <- readsCountsList[[cell1Idx]][!(readsCountsList[[cell1Idx]]$key %in% readsCountsList[[cell0Idx]]$key),]

          # to be in t1, must have data in both cells and have !0 in der alleles
          numT1 <- sum(inBoth[[1]][, derReadsCol] != 0 & inBoth[[2]][, derReadsCol] != 0)
          # to be in t2, must have (data in only cell0 and !0 in der alleles) || (data in both cells and !0 in cell0 der but 0 in cell1 der)
          numT2 <- sum(onlyCell0[, derReadsCol] != 0) + sum(inBoth[[1]][, derReadsCol] != 0 & inBoth[[2]][, derReadsCol] == 0)
          # to be in t3, must have (data in only cell3 and !0 in der alleles) || (data in both cells and !0 in cell1 der but 0 in cell0 der)
          numT3 <- sum(onlyCell1[, derReadsCol] != 0) + sum(inBoth[[1]][, derReadsCol] == 0 & inBoth[[2]][, derReadsCol] != 0)

          data.frame(cell0=cellNames[cell0Idx], cell1=cellNames[cell1Idx], t1=numT1, t2=numT2, t3=numT3)
        }))
      }))
    }))
  }
  rownames(mutCounts) <- NULL

  mutCounts_m <- melt(mutCounts, id.vars=c("cell0", "cell1"))
  mutCounts_m$left <-  sapply(mutCounts_m$cell0, convertCellNumToNewickLabel)
  mutCounts_m$right <- sapply(mutCounts_m$cell1, convertCellNumToNewickLabel)

  merged <- merge(mutCounts_m, treeBranches_m, by=c("left", "right", "variable"))
  colnames(merged) <- c("left", "right", "variable", "cell0", "cell1", "inferred", "true")

  write.table(allIndMutCounts, file=paste0(outputFile, "_indv.txt"), quote=F, sep="\t", col.names=T, row.names=F)
  write.table(merged, file=paste0(outputFile, "_pairs.txt"), quote=F, sep="\t", col.names=T, row.names=F)

  list(indv=allIndMutCounts, paired=merged)
}

saveMutCountsToFile <- function(paramSet, numCells, trueMutCounts, obsMutCounts) {
  trueOutputFile <- paste0(dataDir, paramSet, "/trueMutCounts_c", numCells)
  obsOutputFile <- paste0(dataDir, paramSet, "/obsMutCounts_c", numCells)

  tumorDepths <- read.table(paste0(dataDir, paramSet, "/tumor_depths_", numCells), stringsAsFactors=F)$V1
  cellIDs <- str_extract(tumorDepths, "cancer_cell_[0-9]*.")
  #cellIDs <- gsub(".$", "", cellIDs)

  trueOrderedMutCounts_indv <- NULL
  obsOrderedMutCounts_indv <- NULL
  trueOrderedMutCounts_paired <- NULL
  obsOrderedMutCounts_paired <- NULL
  for(i in 1:(length(tumorDepths)-1)) {
    #cell0_id <- paste0(cellIDs[i], "$")
    cell0_id <- cellIDs[i]
    for(j in (i+1):length(tumorDepths)) {
      #cell1_id <- paste0(cellIDs[j], "$")
      cell1_id <- cellIDs[j]
      currTrueMutCounts <- subset(trueMutCounts[["paired"]], with(trueMutCounts[["paired"]], grepl(gsub(".$", "$", cell0_id), cell0)) & with(trueMutCounts[["paired"]], grepl(gsub(".$", "$", cell1_id), cell1)))[, c("inferred", "variable")]
      currObsMutCounts <- subset(obsMutCounts[["paired"]], with(obsMutCounts[["paired"]], grepl(cell0_id, cell0)) & with(obsMutCounts[["paired"]], grepl(cell1_id, cell1)))[, c("inferred", "variable")]
      trueOrderedMutCounts_paired <- rbind(trueOrderedMutCounts_paired, data.frame(cell0=cell0_id, cell1=cell1_id, trueCount=currTrueMutCounts))
      obsOrderedMutCounts_paired <- rbind(obsOrderedMutCounts_paired, data.frame(cell0=cell0_id, cell1=cell1_id, obsCount=currObsMutCounts))
    }
    currTrueMutCounts <- subset(trueMutCounts[["indv"]], with(trueMutCounts[["indv"]], grepl(gsub(".$", "$", cell0_id), cell)))$t_muts
    currObsMutCounts <- subset(obsMutCounts[["indv"]], with(obsMutCounts[["indv"]], grepl(cell0_id, cell)))$t_muts
    trueOrderedMutCounts_indv <- rbind(trueOrderedMutCounts_indv, data.frame(cell=cell0_id, trueCount=currTrueMutCounts))
    obsOrderedMutCounts_indv <- rbind(obsOrderedMutCounts_indv, data.frame(cell=cell0_id, obsCount=currObsMutCounts))
  }
  write.table(trueOrderedMutCounts_paired[, "trueCount.inferred"], file=paste0(trueOutputFile, "_paired.txt"), quote=F, row.names=F, col.names=F)
  write.table(obsOrderedMutCounts_paired[, "obsCount.inferred"], file=paste0(obsOutputFile, "_paired.txt"), quote=F, row.names=F, col.names=F)
  write.table(paste0("gsl_vector_set(trueMutCounts, ", 0:nrow(trueOrderedMutCounts_paired), ", ", trueOrderedMutCounts_paired[, "trueCount.inferred"], ");"), file=paste0(trueOutputFile, "_paired.cpp"), quote=F, row.names=F, col.names=F)
  write.table(paste0("gsl_vector_set(obsMutCounts, ", 0:nrow(obsOrderedMutCounts_paired), ", ", obsOrderedMutCounts_paired[, "obsCount.inferred"], ");"), file=paste0(obsOutputFile, "_paired.cpp"), quote=F, row.names=F, col.names=F)

  write.table(trueOrderedMutCounts_indv[, "trueCount"], file=paste0(trueOutputFile, "_indv.txt"), quote=F, row.names=F, col.names=F)
  write.table(obsOrderedMutCounts_indv[, "obsCount"], file=paste0(obsOutputFile, "_indv.txt"), quote=F, row.names=F, col.names=F)
  write.table(paste0("gsl_vector_set(trueMutCounts, ", 0:nrow(trueOrderedMutCounts_indv), ", ", trueOrderedMutCounts_indv[, "trueCount"], ");"), file=paste0(trueOutputFile, "_indv.cpp"), quote=F, row.names=F, col.names=F)
  write.table(paste0("gsl_vector_set(obsMutCounts, ", 0:nrow(obsOrderedMutCounts_indv), ", ", obsOrderedMutCounts_indv[, "obsCount"], ");"), file=paste0(obsOutputFile, "_indv.cpp"), quote=F, row.names=F, col.names=F)
}

makeIndvTrueObsInferredCorrPlot <- function(trueIndvMutCounts, obsIndvMutCounts, inferredIndvMutCounts, datasetName) {
  colnames(inferredIndvMutCounts)[colnames(inferredIndvMutCounts) == "true"] <- "branchLength"
  colnames(inferredIndvMutCounts)[colnames(inferredIndvMutCounts) == "value"] <- "numInferred"
  colnames(inferredIndvMutCounts)[colnames(inferredIndvMutCounts) == "variable"] <- "branch"
  colnames(obsIndvMutCounts)[colnames(obsIndvMutCounts) == "true"] <- "branchLength"
  colnames(obsIndvMutCounts)[colnames(obsIndvMutCounts) == "t_muts"] <- "numObserved"
  obsIndvMutCounts$branch <- "t_muts"
  colnames(trueIndvMutCounts)[colnames(trueIndvMutCounts) == "true"] <- "branchLength"
  colnames(trueIndvMutCounts)[colnames(trueIndvMutCounts) == "t_muts"] <- "numTrue"
  trueIndvMutCounts$branch <- "t_muts"

  mergeCols <- c("left", "branch", "branchLength")
  trueAndInferred <- merge(trueIndvMutCounts, inferredIndvMutCounts, by=mergeCols)
  allMutCounts <- merge(trueAndInferred, obsIndvMutCounts, by=mergeCols)

  allMutCounts_m <- melt(allMutCounts[,c(mergeCols, "numTrue", "numObserved", "numInferred")], id.vars=c("left", "branch", "branchLength"))
  lm_eqn <- function(df, formula){
    m <- lm(formula, df)
    rSq <- substitute(italic(r)^2~"="~r2, list(r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(rSq))
  }
  lmStr <- do.call(rbind, lapply(unique(allMutCounts_m$branch), FUN=function(br) {
    do.call(rbind, lapply(unique(allMutCounts_m$variable), FUN=function(var) {
      data.frame(branch=br, variable = var, label=lm_eqn(subset(allMutCounts_m, branch == br & variable == var), value ~ branchLength))
    }))
  }))
  #pNumMutsVsBranches <- ggplot(allMutCounts_m, aes(x=branchLength, y=value, colour=branch)) + geom_point() + facet_grid(variable ~ branch) + theme_bw() + labs(y="num mutations", title="num mutations vs branch length, by observation stage") + geom_text(data=lmStr, aes(x=-Inf, y=Inf, label=label, hjust=0, vjust=1), parse=TRUE, colour="black", show.legend=F) + geom_smooth(method="lm", colour="gray35", show.legend=F)
  pNumMutsVsBranches <- ggplot(allMutCounts_m, aes(x=branchLength, y=value, colour=branch)) + geom_point() + facet_grid(variable ~ branch, scales="free") + theme_bw() + labs(y="num mutations", title=datasetName) + geom_text(data=lmStr, aes(x=-Inf, y=Inf, label=label, hjust=0, vjust=1), parse=TRUE, colour="black", show.legend=F) + geom_smooth(method="lm", colour="gray35", show.legend=F)

  obsVsTrue_m <- melt(allMutCounts[,c(mergeCols, "numTrue", "numObserved")], id.vars=c("left", "branch", "branchLength", "numTrue"))
  infVsTrue_m <- melt(allMutCounts[,c(mergeCols, "numTrue", "numInferred")], id.vars=c("left", "branch", "branchLength", "numTrue"))
  infVsObs_m <- melt(allMutCounts[,c(mergeCols, "numObserved", "numInferred")], id.vars=c("left", "branch", "branchLength", "numObserved"))
  
  obsVsTrue_p <- ggplot(obsVsTrue_m, aes(x=numTrue, y=value, colour=branch)) + labs(y="num observed", title=datasetName) + geom_text(data=ddply(obsVsTrue_m,.(branch), lm_eqn, value ~ numTrue), aes(x=-Inf, y=Inf, label=V1, hjust=0, vjust=1), parse=TRUE, colour="black", show.legend=F)
  infVsTrue_p <- ggplot(infVsTrue_m, aes(x=numTrue, y=value, colour=branch)) + labs(y="num inferred") + geom_text(data=ddply(infVsTrue_m,.(branch), lm_eqn, value ~ numTrue), aes(x=-Inf, y=Inf, label=V1, hjust=0, vjust=1), parse=TRUE, colour="black", show.legend=F)
  infVsObs_p <- ggplot(infVsObs_m, aes(x=numObserved, y=value, colour=branch)) + labs(y="num inferred", x="num observed (by observed read count)") + geom_text(data=ddply(infVsObs_m,.(branch), lm_eqn, value ~ numObserved), aes(x=-Inf, y=Inf, label=V1, hjust=0, vjust=1), parse=TRUE, colour="black", show.legend=F)
  
  plotlist <- list(obsVsTrue_p, infVsTrue_p, infVsObs_p)
  plotlist <- lapply(plotlist, FUN=function(p) {p + geom_point() + facet_wrap(~branch, scales="free") + theme_bw() + geom_smooth(method="lm", colour="gray35", show.legend=F) + theme(legend.position="none")})
  
  pGrid <- plot_grid(plotlist=plotlist, align='vh', nrow=1, byrow=T)

  list(pNumMutsVsBranches=pNumMutsVsBranches, mutCorr=pGrid)
}

makePairedTrueObsInferredCorrPlot <- function(truePairedMutCounts, obsPairedMutCounts, inferredPairedMutCounts, datasetName) {
  colnames(inferredPairedMutCounts)[colnames(inferredPairedMutCounts) == "true"] <- "branchLength"
  colnames(inferredPairedMutCounts)[colnames(inferredPairedMutCounts) == "inferred"] <- "numInferred"
  colnames(inferredPairedMutCounts)[colnames(inferredPairedMutCounts) == "variable"] <- "branch"
  colnames(obsPairedMutCounts)[colnames(obsPairedMutCounts) == "true"] <- "branchLength"
  colnames(obsPairedMutCounts)[colnames(obsPairedMutCounts) == "inferred"] <- "numObserved"
  colnames(obsPairedMutCounts)[colnames(obsPairedMutCounts) == "variable"] <- "branch"
  colnames(truePairedMutCounts)[colnames(truePairedMutCounts) == "true"] <- "branchLength"
  colnames(truePairedMutCounts)[colnames(truePairedMutCounts) == "inferred"] <- "numTrue"
  colnames(truePairedMutCounts)[colnames(truePairedMutCounts) == "variable"] <- "branch"

  mergeCols <- c("left", "right", "branch", "branchLength")
  trueAndInferred <- merge(truePairedMutCounts, inferredPairedMutCounts, by=mergeCols)
  allMutCounts <- merge(trueAndInferred, obsPairedMutCounts, by=mergeCols)

  # get subscripts for labelling
  allMutCounts$branch <- factor(allMutCounts$branch)
  levels(allMutCounts$branch) <- c("t[1]", "t[2]", "t[3]")
  branchColors <- hue_pal()(3) # 3 branches

  allMutCounts_m <- melt(allMutCounts[,c(mergeCols, "numTrue", "numObserved", "numInferred")], id.vars=c("left", "right", "branch", "branchLength"))
  allMutCounts_m$variable <- factor(allMutCounts_m$variable, levels=c("numTrue", "numObserved", "numInferred"))
  levels(allMutCounts_m$variable) <- c("num~true", "num~observed", "num~inferred")
  lm_eqn <- function(df, formula){
    m <- lm(formula, df)
    rSq <- substitute(italic(r)^2~"="~r2, list(r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(rSq))
  }
  lmStr <- do.call(rbind, lapply(levels(allMutCounts_m$branch), FUN=function(br) {
    do.call(rbind, lapply(levels(allMutCounts_m$variable), FUN=function(var) {
      data.frame(branch=br, variable = var, label=lm_eqn(subset(allMutCounts_m, branch == br & variable == var), value ~ branchLength))
    }))
  }))
  lmStr$variable <- factor(lmStr$variable, levels=levels(allMutCounts_m$variable))
  pNumMutsVsBranches <- ggplot(allMutCounts_m, aes(x=branchLength, y=value, colour=branch)) + geom_point(alpha=0.5) + facet_grid(variable ~ branch, scales="free", labeller=label_parsed) + theme_bw() + labs(y="num mutations", x="true branch length", title=datasetName) + geom_text(data=lmStr, aes(x=-Inf, y=Inf, label=label, hjust=0, vjust=1), parse=TRUE, colour="black", show.legend=F) + geom_smooth(method="lm", colour="gray35", show.legend=F) + guides(colour=guide_legend(override.aes=list(alpha=1))) + scale_colour_manual(labels=scales::parse_format(), values=branchColors)

  obsVsTrue_m <- melt(allMutCounts[,c(mergeCols, "numTrue", "numObserved")], id.vars=c("left", "right", "branch", "branchLength", "numTrue"))
  infVsTrue_m <- melt(allMutCounts[,c(mergeCols, "numTrue", "numInferred")], id.vars=c("left", "right", "branch", "branchLength", "numTrue"))
  infVsObs_m <- melt(allMutCounts[,c(mergeCols, "numObserved", "numInferred")], id.vars=c("left", "right", "branch", "branchLength", "numObserved"))
  
  obsVsTrue_p <- ggplot(obsVsTrue_m, aes(x=numTrue, y=value, colour=branch)) + labs(y="num observed", x="num true", title=datasetName) + geom_text(data=ddply(obsVsTrue_m,.(branch), lm_eqn, value ~ numTrue), aes(x=-Inf, y=Inf, label=V1, hjust=0, vjust=1), parse=TRUE, colour="black", show.legend=F)
  infVsTrue_p <- ggplot(infVsTrue_m, aes(x=numTrue, y=value, colour=branch)) + labs(y="num inferred", x="num true") + geom_text(data=ddply(infVsTrue_m,.(branch), lm_eqn, value ~ numTrue), aes(x=-Inf, y=Inf, label=V1, hjust=0, vjust=1), parse=TRUE, colour="black", show.legend=F)
  infVsObs_p <- ggplot(infVsObs_m, aes(x=numObserved, y=value, colour=branch)) + labs(y="num inferred", x="num observed") + geom_text(data=ddply(infVsObs_m,.(branch), lm_eqn, value ~ numObserved), aes(x=-Inf, y=Inf, label=V1, hjust=0, vjust=1), parse=TRUE, colour="black", show.legend=F)
  
  plotlist <- list(obsVsTrue_p, infVsTrue_p, infVsObs_p)
  plotlist <- lapply(plotlist, FUN=function(p) {p + geom_point(alpha=0.5) + facet_wrap(~branch, scales="free", labeller=label_parsed) + theme_bw() + geom_smooth(method="lm", colour="gray35", show.legend=F) + theme(legend.position="none", axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) + guides(colour=guide_legend(override.aes=list(alpha=1))) + scale_colour_manual(labels=scales::parse_format(), values=branchColors)})
  
  pGrid <- plot_grid(plotlist=plotlist, align='vh', nrow=1, byrow=T)

  list(pNumMutsVsBranches=pNumMutsVsBranches, mutCorr=pGrid)
}

makePairedObsInferredCorrPlot <- function(obsPairedMutCounts, inferredPairedMutCounts, datasetName) {
  colnames(inferredPairedMutCounts)[colnames(inferredPairedMutCounts) == "true"] <- "branchLength"
  colnames(inferredPairedMutCounts)[colnames(inferredPairedMutCounts) == "inferred"] <- "numInferred"
  colnames(inferredPairedMutCounts)[colnames(inferredPairedMutCounts) == "variable"] <- "branch"
  colnames(obsPairedMutCounts)[colnames(obsPairedMutCounts) == "true"] <- "branchLength"
  colnames(obsPairedMutCounts)[colnames(obsPairedMutCounts) == "inferred"] <- "numObserved"
  colnames(obsPairedMutCounts)[colnames(obsPairedMutCounts) == "variable"] <- "branch"

  mergeCols <- c("left", "right", "branch", "branchLength")
  allMutCounts <- merge(inferredPairedMutCounts, obsPairedMutCounts, by=mergeCols)

  # get subscripts for labelling
  allMutCounts$branch <- factor(allMutCounts$branch)
  levels(allMutCounts$branch) <- c("t[1]", "t[2]", "t[3]")
  branchColors <- hue_pal()(3) # 3 branches

  lm_eqn <- function(df, formula){
    m <- lm(formula, df)
    rSq <- substitute(italic(r)^2~"="~r2, list(r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(rSq))
  }
  lmStr <- do.call(rbind, lapply(levels(allMutCounts$branch), FUN=function(br) {
      data.frame(branch=br, r2=lm_eqn(subset(allMutCounts, branch == br), numInferred ~ numObserved))
  }))

  infVsObs_m <- melt(allMutCounts[,c(mergeCols, "numObserved", "numInferred")], id.vars=c("left", "right", "branch", "branchLength", "numObserved"))
  
  infVsObs_p <- ggplot(infVsObs_m, aes(x=numObserved, y=value, colour=branch)) + labs(y="num inferred", x="num observed") + geom_text(data=lmStr, aes(x=-Inf, y=Inf, label=r2, hjust=0, vjust=1), parse=TRUE, colour="black", show.legend=F) + geom_point(alpha=0.5) + facet_wrap(~branch, scales="free", labeller=label_parsed) + theme_bw() + geom_smooth(method="lm", colour="gray35", show.legend=F) + theme(legend.position="none") + guides(colour=guide_legend(override.aes=list(alpha=1))) + scale_colour_manual(labels=scales::parse_format(), values=branchColors)
  
  infVsObs_p
}


makePairedTrueObsCorrPlot <- function(truePairedMutCounts, obsPairedMutCounts, datasetName) {
  colnames(obsPairedMutCounts)[colnames(obsPairedMutCounts) == "true"] <- "branchLength"
  colnames(obsPairedMutCounts)[colnames(obsPairedMutCounts) == "inferred"] <- "numObserved"
  colnames(obsPairedMutCounts)[colnames(obsPairedMutCounts) == "variable"] <- "branch"
  colnames(truePairedMutCounts)[colnames(truePairedMutCounts) == "true"] <- "branchLength"
  colnames(truePairedMutCounts)[colnames(truePairedMutCounts) == "inferred"] <- "numTrue"
  colnames(truePairedMutCounts)[colnames(truePairedMutCounts) == "variable"] <- "branch"

  mergeCols <- c("left", "right", "branch", "branchLength")
  allMutCounts <- merge(truePairedMutCounts, obsPairedMutCounts, by=mergeCols)

  allMutCounts_m <- melt(allMutCounts[,c(mergeCols, "numTrue", "numObserved")], id.vars=c("left", "right", "branch", "branchLength"))
  lm_eqn <- function(df, formula){
    m <- lm(formula, df)
    rSq <- substitute(italic(r)^2~"="~r2, list(r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(rSq))
  }
  lmStr <- do.call(rbind, lapply(unique(allMutCounts_m$branch), FUN=function(br) {
    do.call(rbind, lapply(unique(allMutCounts_m$variable), FUN=function(var) {
      data.frame(branch=br, variable = var, label=lm_eqn(subset(allMutCounts_m, branch == br & variable == var), value ~ branchLength))
    }))
  }))
  pNumMutsVsBranches <- ggplot(allMutCounts_m, aes(x=branchLength, y=value, colour=branch)) + geom_point(alpha=.2) + facet_grid(variable ~ branch, scales="free") + theme_bw() + labs(y="num mutations", title=datasetName) + geom_text(data=lmStr, aes(x=-Inf, y=Inf, label=label, hjust=0, vjust=1), parse=TRUE, colour="black", show.legend=F) + geom_smooth(method="lm", colour="gray35", show.legend=F)

  obsVsTrue_m <- melt(allMutCounts[,c(mergeCols, "numTrue", "numObserved")], id.vars=c("left", "right", "branch", "branchLength", "numTrue"))
  
  obsVsTrue_p <- ggplot(obsVsTrue_m, aes(x=numTrue, y=value, colour=branch)) + labs(y="num observed", title=datasetName) + geom_text(data=ddply(obsVsTrue_m,.(branch), lm_eqn, value ~ numTrue), aes(x=-Inf, y=Inf, label=V1, hjust=0, vjust=1), parse=TRUE, colour="black", show.legend=F)
  
  plotlist <- list(obsVsTrue_p)
  plotlist <- lapply(plotlist, FUN=function(p) {p + geom_point() + facet_wrap(~branch, scales="free") + theme_bw() + geom_smooth(method="lm", colour="gray35", show.legend=F) + theme(legend.position="none")})
  pGrid <- plot_grid(plotlist=plotlist, align='vh', nrow=1, byrow=T)

  list(pNumMutsVsBranches=pNumMutsVsBranches, mutCorr=pGrid)
}




# based on counts from running inference
getInferredPerBranchMutCounts <- function(newickString, branchLength, paramSet, forceRecalc=F) {
  outputFile <- paste0(hmmFile, ".inferredMutCounts")
  if(!forceRecalc && file.exists(paste0(outputFile, "_indv.txt")) && file.exists(paste0(outputFile, "_pairs.txt"))) {
    indvDat <- read.table(paste0(outputFile, "_indv.txt"), sep="\t", header=T)
    pairDat <- read.table(paste0(outputFile, "_pairs.txt"), sep="\t", header=T)
    return(list(indv=indvDat, paired=pairDat))
  }

  treeObj <- read.tree(text=newickString)
  mrcaMat <- mrca(treeObj)
  distMat <- dist.nodes(treeObj)
  nodeDepths <- node.depth.edgelength(treeObj)

  # left | right | t1 | t2 | t3
  # left/right correspond to node labels in the newick file
  treeBranches <- calcTreeBranches(mrcaMat, distMat, nodeDepths, branchLength)
  treeBranches_m <- melt(treeBranches, id.vars=c("left", "right"))

  # first check if this run finished correctly, and the string FINAL HMM appears in hmmFile
  if(suppressWarnings(length(system(paste0("/bin/bash -c ", shQuote(sprintf("grep \"FINAL HMM\" %s", hmmFile))), intern=T)) == 0)) {
    return(NULL)
  }

  # allIndEstMutCounts
  indHmmNames <- system(paste0("/bin/bash -c ", shQuote(sprintf("grep \"^HMM\" %s | grep -v \",\" | uniq | sed -e 's/HMM 0\\s\\+//' -e 's/(//' -e 's/)//' -e 's/://' ", hmmFile))), intern=T)
  allIndEstMutCounts <- read.table(text=system(paste0("/bin/bash -c ", shQuote(sprintf("sed '1,/FINAL HMM/d' %s | tac | sed '/allIndEstMutCounts/q' | tac | sed '/allPairedEstMutCounts/q' | tail -n +2 | head -n -2", hmmFile))), intern=T))
  colnames(allIndEstMutCounts) <- "value"
  allIndEstMutCounts$cell <- indHmmNames
  allIndEstMutCounts$variable <- "t_muts"
  allIndEstMutCounts$left <- sapply(allIndEstMutCounts$cell, convertCellNumToNewickLabel)
  #allIndEstMutCounts$true <- nodeDepths[allIndEstMutCounts$left] + branchLength
  allIndEstMutCounts$true <- nodeDepths[allIndEstMutCounts$left]# + branchLength Thu 30 Jun 2022 03:01:45 PM PDT I think this might be wrong?

  # cell0 | cell1 | t1
  # cell0 | cell1 | t2
  # cell0 | cell1 | t3
  # ...
  # need to separate cell name reading from param reading in case doing nearest n
  hmmLines <- system(paste0("/bin/bash -c ", shQuote(sprintf("grep \"^HMM\" %s | grep \",\" | sed -e 's/HMM\\s\\+//' -e 's/(//' -e 's/)//' -e 's/,/\t/' -e 's/://' -e 's/\\s\\+/\t/'", hmmFile))), intern=T)
  if(length(hmmLines) == 0) {
    return(NULL)
  }
  analyzedHmmNames <- read.table(text=hmmLines)
  colnames(analyzedHmmNames) <- c("hmmIdx", "cell0", "cell1")
  paramIdxToExtract <- do.call(c, lapply(analyzedHmmNames$hmmIdx, FUN=function(i) {c(3*i, (3*i)+1, (3*i)+2)})) # don't extract hmms we skipped
  paramIdxToExtract <- paramIdxToExtract + 1 # paramsToEst is 0 indexed, R vector is 1 indexed
  #unlabeledHmmParamsRaw <- read.table(text=system(paste0("/bin/bash -c ", shQuote(sprintf("sed '1,/FINAL HMM/d' %s | tac | sed '/paramsToEst/q' | tac | sed '/fixedParams/q' | tail -n +2 | head -n -2", hmmFile))), intern=T))
  #hmmParamsRaw <- data.frame(cell0=do.call(c, lapply(analyzedHmmNames$cell0, rep, 3)), cell1=do.call(c, lapply(analyzedHmmNames$cell1, rep, 3)), value=unlabeledHmmParamsRaw[paramIdxToExtract,])
  #hmmParamsRaw$variable <- c("t1_muts", "t2_muts", "t3_muts")

  # allPairedEstMutCounts
  unlabeledAllPairedEstMutCounts <- read.table(text=system(paste0("/bin/bash -c ", shQuote(sprintf("sed '1,/FINAL HMM/d' %s | tac | sed '/allPairedEstMutCounts/q' | tac | sed '/IndThenPairs2Stages3TrParam2DegPolyHMMWithMuts/q' | tail -n +2 | head -n -2", hmmFile))), intern=T))
  allPairedEstMutCounts <- data.frame(cell0=do.call(c, lapply(analyzedHmmNames$cell0, rep, 3)), cell1=do.call(c, lapply(analyzedHmmNames$cell1, rep, 3)), value=unlabeledAllPairedEstMutCounts[paramIdxToExtract,])
  allPairedEstMutCounts$variable <- c("t1", "t2", "t3")

  allPairedEstMutCounts$left <-  sapply(allPairedEstMutCounts$cell0, convertCellNumToNewickLabel)
  allPairedEstMutCounts$right <- sapply(allPairedEstMutCounts$cell1, convertCellNumToNewickLabel)
  
  allPairedEstMutCountsMerged <- merge(allPairedEstMutCounts, treeBranches_m, by=c("left", "right", "variable"))
  colnames(allPairedEstMutCountsMerged) <- c("left", "right", "variable", "cell0", "cell1", "inferred", "true")

  write.table(allIndEstMutCounts, file=paste0(outputFile, "_indv.txt"), quote=F, sep="\t", col.names=T, row.names=F)
  write.table(allPairedEstMutCountsMerged, file=paste0(outputFile, "_pairs.txt"), quote=F, sep="\t", col.names=T, row.names=F)
  list(indv=allIndEstMutCounts, paired=allPairedEstMutCountsMerged)
}


# not per branch
compareTrueAndObsMutCounts <- function(paramSet, forceRecalc=F) {
  outputFile <- paste0(dataDir, paramSet, "/trueVsObsMutCounts.txt")
  if(!forceRecalc && file.exists(outputFile)) {
    merged <- read.table(outputFile, sep="\t", header=T, stringsAsFactors=F)
    return(merged)
  }
  snpSummariesFileList <- system(paste0("find ", dataDir, paramSet, " -maxdepth 1 -name \"snpSummaryCN_cancer_cell_*\" -not -name \"*err\"| sort -V"), intern=T)
  snpSummariesList <- lapply(snpSummariesFileList, FUN=function(f) {
    cellName <- basename(f)
    tab <- read.table(f, sep="\t", header=F)
    colnames(tab) <- c("binIdx", "snpIdx", "ancAlleles", "derAlleles", "ancReads", "derReads", "CN")
    tab$key <- paste0(tab$chr, ":", tab$start, "-", tab$end)
    tab
  })
  cellNames <- sapply(snpSummariesFileList, basename)
  names(snpSummariesList) <- cellNames

  # for each cell
  mutTypeCounts <- do.call(rbind, lapply(1:(length(cellNames)), FUN=function(cellIdx) {
    # remove sites with 0 reads
    allSites <- snpSummariesList[[cellIdx]]
    sitesWithReads <- allSites[with(allSites, ancReads > 0 | derReads > 0),]

    # expected ancestral homozygous: only ancestral alleles
    expAncHom <- sum(with(sitesWithReads, ancAlleles > 0 & derAlleles == 0))
    expHet <- sum(with(sitesWithReads, ancAlleles > 0 & derAlleles > 0))
    expDerHom <- sum(with(sitesWithReads, ancAlleles == 0 & derAlleles > 0))

    obsAncHom <- sum(with(sitesWithReads, ancReads > 0 & derReads == 0))
    obsHet <- sum(with(sitesWithReads, ancReads > 0 & derReads > 0))
    obsDerHom <- sum(with(sitesWithReads, ancReads == 0 & derReads > 0))
    data.frame(cell=cellNames[cellIdx], expAncestralHom=expAncHom, expHet=expHet, expDerivedHom=expDerHom, obsAncestralHom=obsAncHom, obsHet=obsHet, obsDerivedHom=obsDerHom)
  }))
  rownames(mutTypeCounts) <- NULL

  write.table(mutTypeCounts, file=outputFile, quote=F, sep="\t", col.names=T, row.names=F)
  mutTypeCounts
}
makeTrueVsObsBoxplots <- function(mutTypeCounts, plotTitle) {
  mutTypeCounts_m <- melt(mutTypeCounts, id.vars="cell")
  mutTypeCounts_m$source <- as.factor(substr(mutTypeCounts_m$variable, 1, 3))
  mutTypeCounts_m$genotype <- as.factor(substr(mutTypeCounts_m$variable, 4, 30))
  p <- ggplot(mutTypeCounts_m, aes(x=source, y=value, colour=genotype)) + geom_boxplot() + facet_wrap(~genotype, scales="free") + theme_bw() #+ stat_compare_means(aes(group=genotype), comparisons=list(c("exp", "obs")), paired=T, method="t.test")
  if(!is.na(plotTitle)) {
    p <- p + labs(title=plotTitle)
  }
  p
}

getTruePerBranchCNAcounts <- function(newickString, branchLength, paramSet, forceRecalc=F, numCellsList=NULL) {
  numCellsFilt <- ifelse(is.null(numCellsList), "", paste0("_c", paste0(numCellsList, collapse="-c")))
  outputFile <- paste0(dataDir, paramSet, "/truePerBranchCNAcounts", numCellsFilt, ".txt")
  if(!forceRecalc && file.exists(outputFile)) {
    merged <- read.table(outputFile, sep="\t", header=T, stringsAsFactors=F)
    return(merged)
  }
  if(is.null(numCellsList)) {
    tumorListFile <- paste0(dataDir, paramSet, "/tumor_depths")
    tumorList <- read.table(tumorListFile)
    tumorTrueBeds <- sapply(tumorList, FUN=function(x) {gsub("depth", "bed", gsub("simu", "true", x))})
    if(!file.exists(tumorTrueBeds[1])) {
      tumorTrueBeds <- sapply(tumorList, FUN=function(x) {gsub(".depth", ".bed", gsub("readDepth", "copyNumber", gsub("simu", "true", x)))})
    }
     
    allDepths <- do.call(cbind, lapply(tumorTrueBeds, FUN=function(x) {
      bed <- read.table(x, sep="\t", header=F, stringsAsFactors=F)
      colnames(bed) <- c("chr", "start", "end", basename(x))
      bed[,4, drop=F]
    }))
    cellNames <- colnames(allDepths)
    
    # use RLE length as an approximation of boundaries of copy number events
    # names are already newick labs
    sharedBreakpoints <- do.call(rbind, lapply(1:(length(cellNames)-1), FUN=function(cell0Idx) {
      do.call(rbind, lapply((cell0Idx+1):length(cellNames), FUN=function(cell1Idx) {
        # use indices of breakpoints as an approximation of shared CNAs
        cell0Breakpoints <- cumsum(rle(allDepths[,cell0Idx])$lengths)
        cell1Breakpoints <- cumsum(rle(allDepths[,cell1Idx])$lengths)
        # num shared breakpoints (t1)
        numT1 <- sum(cell0Breakpoints %in% cell1Breakpoints)
        # num breakpoints in only cell0 (t2)
        numT2 <- sum(!(cell0Breakpoints %in% cell1Breakpoints))
        # num breakpoints in only cell1 (t3)
        numT3 <- sum(!(cell1Breakpoints %in% cell0Breakpoints))
        data.frame(cell0=cellNames[cell0Idx], cell1=cellNames[cell1Idx], t1=numT1, t2=numT2, t3=numT3)
      }))
    }))
  } else {
    sharedBreakpoints <- do.call(rbind, lapply(numCellsList, FUN=function(numCells) {
      tumorListFile <- paste0(dataDir, paramSet, "/tumor_depths_", numCells)
      tumorList <- read.table(tumorListFile)
      tumorTrueBeds <- sapply(tumorList, FUN=function(x) {gsub("depth", "bed", gsub("simu", "true", x))})
      if(!file.exists(tumorTrueBeds[1])) {
        tumorTrueBeds <- sapply(tumorList, FUN=function(x) {gsub(".depth", ".bed", gsub("readDepth", "copyNumber", gsub("simu", "true", x)))})
      }
       
      allDepths <- do.call(cbind, lapply(tumorTrueBeds, FUN=function(x) {
        bed <- read.table(x, sep="\t", header=F, stringsAsFactors=F)
        colnames(bed) <- c("chr", "start", "end", basename(x))
        bed[,4, drop=F]
      }))
      cellNames <- colnames(allDepths)
      
      # use RLE length as an approximation of boundaries of copy number events
      # names are already newick labs
      do.call(rbind, lapply(1:(length(cellNames)-1), FUN=function(cell0Idx) {
        do.call(rbind, lapply((cell0Idx+1):length(cellNames), FUN=function(cell1Idx) {
          # use indices of breakpoints as an approximation of shared CNAs
          cell0Breakpoints <- cumsum(rle(allDepths[,cell0Idx])$lengths)
          cell1Breakpoints <- cumsum(rle(allDepths[,cell1Idx])$lengths)
          # num shared breakpoints (t1)
          numT1 <- sum(cell0Breakpoints %in% cell1Breakpoints)
          # num breakpoints in only cell0 (t2)
          numT2 <- sum(!(cell0Breakpoints %in% cell1Breakpoints))
          # num breakpoints in only cell1 (t3)
          numT3 <- sum(!(cell1Breakpoints %in% cell0Breakpoints))
          data.frame(cell0=cellNames[cell0Idx], cell1=cellNames[cell1Idx], t1=numT1, t2=numT2, t3=numT3)
        }))
      }))
    }))
  }

  sharedBreakpoints_m <- melt(sharedBreakpoints, id.vars=c("cell0", "cell1"))
  sharedBreakpoints_m$left <-  sapply(sharedBreakpoints_m$cell0, convertCellNumToNewickLabel)
  sharedBreakpoints_m$right <- sapply(sharedBreakpoints_m$cell1, convertCellNumToNewickLabel)

  # left | right | t1 | t2 | t3
  # left/right correspond to node labels in the newick file
  treeObj <- read.tree(text=newickString)
  mrcaMat <- mrca(treeObj)
  distMat <- dist.nodes(treeObj)
  nodeDepths <- node.depth.edgelength(treeObj)
  treeBranches <- calcTreeBranches(mrcaMat, distMat, nodeDepths, branchLength)
  treeBranches_m <- melt(treeBranches, id.vars=c("left", "right"))

  merged <- merge(sharedBreakpoints_m, treeBranches_m, by=c("left", "right", "variable"))
  colnames(merged) <- c("left", "right", "variable", "cell0", "cell1", "inferred", "true")

  write.table(merged, file=outputFile, quote=F, sep="\t", col.names=T, row.names=F)
  merged
}

# returns matrix of true copy numbers from true_cancer_cell*, rows are cells, cols are genomic positions
readAllTumorTrueBeds <- function(dataDir, paramSet, numCells) {
  tumorListFile <- paste0(dataDir, paramSet, "/tumor_depths")
  if(!is.null(numCells)) {
    tumorListFile <- paste0(tumorListFile, "_", numCells)
  }
  
  if(!file.exists(tumorListFile)) {
    next
  }
  tumorList <- read.table(tumorListFile)
  tumorTrueBeds <- sapply(tumorList, FUN=function(x) {gsub("depth", "bed", gsub("simu", "true", x))})
  if(!file.exists(tumorTrueBeds[1])) {
    tumorTrueBeds <- sapply(tumorList, FUN=function(x) {gsub(".depth", ".bed", gsub("readDepth", "copyNumber", gsub("simu", "true", x)))})
  }
  
  tumorDepths <- lapply(tumorTrueBeds, FUN=function(x) {
    bed <- read.table(x, sep="\t", header=F, stringsAsFactors=F)
    colnames(bed) <- c("chr", "start", "end", "cn")
    bed$cellFile <- str_extract(x, "cancer_cell_[0-9]*.")
    bed$newickLab <- factor(128 - as.numeric(gsub(".$", "", gsub("cancer_cell_", "", bed$cellFile))))
    #bed$newickLab <- factor(1 + as.numeric(gsub(".$", "", gsub("cancer_cell_", "", bed$cellFile))))
    bed
  })
  
  newickLabs <- factor(128 - as.numeric(gsub(".$", "", gsub("cancer_cell_", "",str_extract(tumorList$V1, "cancer_cell_[0-9]*.")))))
  #newickLabs <- factor(1 + as.numeric(gsub(".$", "", gsub("cancer_cell_", "",str_extract(tumorList$V1, "cancer_cell_[0-9]*."))))) # decided cells based on cell file numbers, for cextremes sets it's easier to use these numbers
  
  allDepths <- do.call(rbind, lapply(tumorDepths, FUN=function(x) {x$cn}))
  rownames(allDepths) <- newickLabs

  allDepths
}

# returns matrix of inferred copy numbers from outputType ("__sconce__*.bed", "*__mean.bed", "*__median.bed", "*__mode.bed"), rows are cells, cols are genomic positions
readAllOutputBedFiles <- function(dataDir, paramSet, key, numCells, outputType) {
  outputFilePrefix <- paste0("output_", key, "_", gsub("/", "_", paramSet), "_k", k, "_c", numCells) # based on scAllP_*sh outBase variable
  if(grepl("Mut", key) & grepl("params[23]", paramSet)) {
    #outputFilePrefix <- paste0("output_", key, "_reuseMutEsts_shortcut_", gsub("/", "_", paramSet), "_k", k, "_c", numCells) # based on scAllP_*sh outBase variable
  }
  unfiltFileList <- system(paste0("find ", dataDir, paramSet, " -maxdepth 1 -name \"", outputFilePrefix, outputType, "\" | sort -V"), intern=T)
  tumorDepths <- read.table(paste0(dataDir, paramSet, "/tumor_depths_", numCells), stringsAsFactors=F)$V1
  cellIDs <- str_extract(tumorDepths, "cancer_cell_[0-9]*.")
  filtFileList <- sapply(cellIDs, FUN=function(cellID) {unfiltFileList[grepl(cellID, unfiltFileList, fixed=T)]})
  dat <- do.call(rbind, lapply(cellIDs, FUN=function(cellID) {
    bed <- read.table(filtFileList[grepl(cellID, filtFileList, fixed=T)], stringsAsFactors=F, sep="\t", header=F)
    bed$V4
  }))
  newickLabs <- factor(128 - as.numeric(gsub(".$", "", gsub("cancer_cell_", "",str_extract(tumorDepths, "cancer_cell_[0-9]*.")))))
  rownames(dat) <- newickLabs
  dat
}

readcnp2cnpFile <- function(cnp2cnpFile) {
  dat <- read.table(cnp2cnpFile, skip=1, row.names=1, stringsAsFactors=F)
  colnames(dat) <- rownames(dat) <- as.character(sapply(rownames(dat), convertCellNumToNewickLabel))
  as.matrix(dat)
}

###########################################
# Distance matrix functions
###########################################
# given mergedTreeBranches, creates a distance matrix based on t2+t3 values for uniqNodes
createInferred_t2_t3_distMat <- function(mergedTreeBranches, uniqNodes) {
  #inferredDistMat <- matrix(rep(0, length(uniqNodes)^2), ncol=length(uniqNodes)) # values from .hmm file
  inferredDistMat <- matrix(rep(NA, length(uniqNodes)^2), ncol=length(uniqNodes)) # values from .hmm file
  diag(inferredDistMat) <- 0
  dimnames(inferredDistMat) <- list(uniqNodes, uniqNodes)
  t2_t3 <- subset(mergedTreeBranches, variable == "t2_t3") # t2+t3 values from .hmm file
  for(i in 1:nrow(t2_t3)) {
    left <- as.character(t2_t3[i,"left"])
    right <- as.character(t2_t3[i,"right"])
    inferredDistMat[left, right] <- t2_t3[t2_t3$left == left & t2_t3$right == right, "inferred"]
    inferredDistMat[right, left] <- t2_t3[t2_t3$left == left & t2_t3$right == right, "inferred"]
  }
  inferredDistMat
}
createInferred_t1_distMat <- function(mergedTreeBranches, uniqNodes) {
  #inferredDistMat <- matrix(rep(0, length(uniqNodes)^2), ncol=length(uniqNodes)) # values from .hmm file
  inferredDistMat <- matrix(rep(NA, length(uniqNodes)^2), ncol=length(uniqNodes)) # values from .hmm file
  diag(inferredDistMat) <- 0
  dimnames(inferredDistMat) <- list(uniqNodes, uniqNodes)
  t1 <- subset(mergedTreeBranches, variable == "t1") # t1 values from .hmm file
  for(i in 1:nrow(t1)) {
    left <- as.character(t1[i,"left"])
    right <- as.character(t1[i,"right"])
    # convert similarity to distance (t1 is shared history)
    inferredDistMat[left, right] <- 1 - (t1[t1$left == left & t1$right == right, "inferred"] / max(t1$inferred))
    inferredDistMat[right, left] <- 1 - (t1[t1$left == left & t1$right == right, "inferred"] / max(t1$inferred))
  }
  inferredDistMat
}

# given mergedTreeBranches, creates a distance matrix based on true values (from the newick file) for uniqNodes
createTrueDistMat <- function(mergedTreeBranches, uniqNodes) {
  trueDistMat <- matrix(rep(0, length(uniqNodes)^2), ncol=length(uniqNodes)) # values from newick file
  dimnames(trueDistMat) <- list(uniqNodes, uniqNodes)
  t2_t3 <- subset(mergedTreeBranches, variable == "t2_t3") # t2+t3 values from .hmm file
  
  for(i in 1:nrow(t2_t3)) {
    left <- as.character(t2_t3[i,"left"])
    right <- as.character(t2_t3[i,"right"])
    trueDistMat[left, right] <- t2_t3[t2_t3$left == left & t2_t3$right == right, "true"]
    trueDistMat[right, left] <- t2_t3[t2_t3$left == left & t2_t3$right == right, "true"]
  }
  trueDistMat
}

# given a cnp2cnp output file name, reads it in, converts it to a dist mat
createCnpDistMat <- function(cnpFilename, cnpRevFilename) {
  if(!file.exists(cnpFilename)) {
    warning(paste0("could not find cnp2cnp files", cnpFilename, ", ", cnpRevFilename, ". Did you run scripts/cnaBedToFasta_wrapper.sh?"))
    return(NULL)
  }
  cnpMat <- readcnp2cnpFile(cnpFilename)
  cnpRevMat <- readcnp2cnpFile(cnpRevFilename)
  cnpMat <- cnpMat + cnpRevMat[rownames(cnpMat), colnames(cnpMat)]
  cnpMat
}

###########################################
# Neighbor joining and tree calculation functions
###########################################
createNJtrees <- function(mergedTreeBranches, paramSet, numCells, key) {
  uniqNodes <- as.character(sort(unique(c(mergedTreeBranches$left, mergedTreeBranches$right))))
  inferredDistMat <- createInferred_t2_t3_distMat(mergedTreeBranches, uniqNodes)
  trueDistMat <- createTrueDistMat(mergedTreeBranches, uniqNodes)

  # get cnp2cnp values
  trueCnpBase <- paste0(dataDir, "/", paramSet, "/tumor_depths_", numCells)
  cnpTrueMat <- createCnpDistMat(paste0(trueCnpBase, "_true.cnp2cnp"), paste0(trueCnpBase, "_trueRev.cnp2cnp"))

  cnpBase <- paste0(dataDir, "/", paramSet, "/output_", key, "_", gsub("/", "_", paramSet), "_k", k, "_c", numCells)
  cnpSconceMat <- createCnpDistMat(paste0(cnpBase, "_roundedsconce.cnp2cnp"), paste0(cnpBase, "_roundedsconceRev.cnp2cnp"))
  if(is.null(cnpSconceMat)) {
    return(NULL)
  }
  cnpMeanMat <- createCnpDistMat(paste0(cnpBase, "_roundedmean.cnp2cnp"), paste0(cnpBase, "_roundedmeanRev.cnp2cnp"))
  cnpMedianMat <- createCnpDistMat(paste0(cnpBase, "_roundedmedian.cnp2cnp"), paste0(cnpBase, "_roundedmedianRev.cnp2cnp"))
  cnpModeMat <- createCnpDistMat(paste0(cnpBase, "_roundedmode.cnp2cnp"), paste0(cnpBase, "_roundedmodeRev.cnp2cnp"))

  cnpTrueTree   <- nj(cnpTrueMat) # tree built on cnp2cnp metric btn true copy number profiles
  cnpSconceTree <- nj(cnpSconceMat) # tree built on cnp2cnp metric btn sconce copy number profiles
  cnpMeanTree   <- nj(cnpMeanMat) # tree built on cnp2cnp metric btn mean copy number profiles
  cnpMedianTree <- nj(cnpMedianMat) # tree built on cnp2cnp metric btn median copy number profiles
  cnpModeTree   <- nj(cnpModeMat) # tree built on cnp2cnp metric btn mode copy number profiles

  # zzs metric
  trueZzsBase <- paste0(dataDir, "/", paramSet, "/tumor_depths_", numCells)
  zzsTrueMat <- createCnpDistMat(paste0(trueZzsBase, "_true.zzs.cnp2cnp"), paste0(trueZzsBase, "_trueRev.zzs.cnp2cnp"))
  zzsBase <- paste0(dataDir, "/", paramSet, "/output_", key, "_", gsub("/", "_", paramSet), "_k", k, "_c", numCells)
  zzsSconceMat <- createCnpDistMat(paste0(zzsBase, "_roundedsconce.zzs.cnp2cnp"), paste0(zzsBase, "_roundedsconceRev.zzs.cnp2cnp"))
  zzsMeanMat <- createCnpDistMat(paste0(zzsBase, "_roundedmean.zzs.cnp2cnp"), paste0(zzsBase, "_roundedmeanRev.zzs.cnp2cnp"))
  zzsMedianMat <- createCnpDistMat(paste0(zzsBase, "_roundedmedian.zzs.cnp2cnp"), paste0(zzsBase, "_roundedmedianRev.zzs.cnp2cnp"))
  zzsModeMat <- createCnpDistMat(paste0(zzsBase, "_roundedmode.zzs.cnp2cnp"), paste0(zzsBase, "_roundedmodeRev.zzs.cnp2cnp"))

  zzsTrueTree   <- nj(zzsTrueMat) # tree built on zzs cnp2cnp metric btn true copy number profiles
  zzsSconceTree <- nj(zzsSconceMat) # tree built on zzs cnp2cnp metric btn sconce copy number profiles
  zzsMeanTree   <- nj(zzsMeanMat) # tree built on zzs cnp2cnp metric btn mean copy number profiles
  zzsMedianTree <- nj(zzsMedianMat) # tree built on zzs cnp2cnp metric btn median copy number profiles
  zzsModeTree   <- nj(zzsModeMat) # tree built on zzs cnp2cnp metric btn mode copy number profiles

  fullTrueNewickTree <- read.tree(text=newickStrings[paramSet]) # tree read directly from newick file
  trueDistTree <- keep.tip(fullTrueNewickTree, uniqNodes) # tree built directly from newick file, filtered for tips in this subset
  inferredDistTree <- njs(inferredDistMat) # tree built on inferred t2+t3 branch lengths, using njs for missing values if nearest10

  # manually calculate euclidean distance
  # get output from our bed files
  trueCNs <- readAllTumorTrueBeds(dataDir, paramSet, numCells) # values from true_cancer_cell* files
  sconceCNs <- readAllOutputBedFiles(dataDir, paramSet, key, numCells, "*__sconce__*.bed") # values from output*bed files
  meanCNs <- readAllOutputBedFiles(dataDir, paramSet, key, numCells, "*__mean.bed")
  medianCNs <- readAllOutputBedFiles(dataDir, paramSet, key, numCells, "*__median.bed")
  modeCNs <- readAllOutputBedFiles(dataDir, paramSet, key, numCells, "*__mode.bed")
  eucTrueTree   <- nj(dist(trueCNs, method="euclidean")) # tree built on euclidean distance btn true copy number profiles
  eucSconceTree <- nj(dist(sconceCNs, method="euclidean")) # tree built on euclidean distance btn sconce copy number profiles
  eucMeanTree   <- nj(dist(meanCNs, method="euclidean")) # tree built on euclidean distance btn mean copy number profiles
  eucMedianTree <- nj(dist(medianCNs, method="euclidean")) # tree built on euclidean distance btn median copy number profiles
  eucModeTree   <- nj(dist(modeCNs, method="euclidean")) # tree built on euclidean distance btn mode copy number profiles

  #list(eucTrueTree=eucTrueTree, eucSconceTree=eucSconceTree, eucMeanTree=eucMeanTree, eucMedianTree=eucMedianTree, eucModeTree=eucModeTree, cnpTrueTree=cnpTrueTree, cnpSconceTree=cnpSconceTree, cnpMeanTree=cnpMeanTree, cnpMedianTree=cnpMedianTree, cnpModeTree=cnpModeTree, trueDistTree=trueDistTree, inferredDistTree=inferredDistTree)
  list(eucTrueTree=eucTrueTree, eucSconceTree=eucSconceTree, eucMeanTree=eucMeanTree, eucMedianTree=eucMedianTree, eucModeTree=eucModeTree, cnpTrueTree=cnpTrueTree, cnpSconceTree=cnpSconceTree, cnpMeanTree=cnpMeanTree, cnpMedianTree=cnpMedianTree, cnpModeTree=cnpModeTree, zzsTrueTree=zzsTrueTree, zzsSconceTree=zzsSconceTree, zzsMeanTree=zzsMeanTree, zzsMedianTree=zzsMedianTree, zzsModeTree=zzsModeTree, trueDistTree=trueDistTree, inferredDistTree=inferredDistTree)
}

createNJtreesCompareMuts <- function(mergedTreeBranches, paramSet, numCells, key) {
  uniqNodes <- as.character(sort(unique(c(mergedTreeBranches$left, mergedTreeBranches$right))))
  inferredDistMat_t2_t3 <- createInferred_t2_t3_distMat(mergedTreeBranches, uniqNodes)
  inferredDistMat_t1 <- createInferred_t1_distMat(mergedTreeBranches, uniqNodes)
  trueDistMat <- createTrueDistMat(mergedTreeBranches, uniqNodes)

  # get cnp2cnp values
  trueCnpBase <- paste0(dataDir, "/", paramSet, "/tumor_depths_", numCells)
  cnpTrueMat <- createCnpDistMat(paste0(trueCnpBase, "_true.cnp2cnp"), paste0(trueCnpBase, "_trueRev.cnp2cnp"))

  cnpBase <- paste0(dataDir, "/", paramSet, "/output_", key, "_", gsub("/", "_", paramSet), "_k", k, "_c", numCells)
  if(grepl("Mut", key) & grepl("params[23]", paramSet)) {
    #cnpBase <- paste0(dataDir, "/", paramSet, "/output_", key, "_reuseMutEsts_shortcut_", gsub("/", "_", paramSet), "_k", k, "_c", numCells)
  }
  cnpSconceMat <- createCnpDistMat(paste0(cnpBase, "_roundedsconce.cnp2cnp"), paste0(cnpBase, "_roundedsconceRev.cnp2cnp"))
  cnpMeanMat <- createCnpDistMat(paste0(cnpBase, "_roundedmean.cnp2cnp"), paste0(cnpBase, "_roundedmeanRev.cnp2cnp"))
  if(is.null(cnpSconceMat)) {
    return(NULL)
  }
  cnpTrueTree   <- nj(cnpTrueMat) # tree built on cnp2cnp metric btn true copy number profiles
  cnpSconceTree <- nj(cnpSconceMat) # tree built on cnp2cnp metric btn sconce copy number profiles
  cnpMeanTree   <- nj(cnpMeanMat) # tree built on cnp2cnp metric btn mean copy number profiles

  # zzs metric
  trueZzsBase <- paste0(dataDir, "/", paramSet, "/tumor_depths_", numCells)
  zzsTrueMat <- createCnpDistMat(paste0(trueZzsBase, "_true.zzs.cnp2cnp"), paste0(trueZzsBase, "_trueRev.zzs.cnp2cnp"))
  zzsBase <- paste0(dataDir, "/", paramSet, "/output_", key, "_", gsub("/", "_", paramSet), "_k", k, "_c", numCells)
  if(grepl("Mut", key) & grepl("params[23]", paramSet)) {
    #zzsBase <- paste0(dataDir, "/", paramSet, "/output_", key, "_reuseMutEsts_shortcut_", gsub("/", "_", paramSet), "_k", k, "_c", numCells)
  }
  zzsSconceMat <- createCnpDistMat(paste0(zzsBase, "_roundedsconce.zzs.cnp2cnp"), paste0(zzsBase, "_roundedsconceRev.zzs.cnp2cnp"))
  zzsMeanMat <- createCnpDistMat(paste0(zzsBase, "_roundedmean.zzs.cnp2cnp"), paste0(zzsBase, "_roundedmeanRev.zzs.cnp2cnp"))

  zzsTrueTree   <- nj(zzsTrueMat) # tree built on zzs cnp2cnp metric btn true copy number profiles
  zzsSconceTree <- nj(zzsSconceMat) # tree built on zzs cnp2cnp metric btn sconce copy number profiles
  zzsMeanTree   <- nj(zzsMeanMat) # tree built on zzs cnp2cnp metric btn mean copy number profiles

  fullTrueNewickTree <- read.tree(text=newickStrings[paramSet]) # tree read directly from newick file
  trueDistTree <- keep.tip(fullTrueNewickTree, uniqNodes) # tree built directly from newick file, filtered for tips in this subset
  inferredDistTree_t2_t3 <- nj(inferredDistMat_t2_t3) # tree built on inferred t2+t3 branch lengths, using njs for missing values if nearest10
  inferredDistTree_t1 <- nj(inferredDistMat_t1) # tree built on inferred t1 branch lengths, using njs for missing values if nearest10

  # manually calculate euclidean distance
  # get output from our bed files
  trueCNs <- readAllTumorTrueBeds(dataDir, paramSet, numCells) # values from true_cancer_cell* files
  sconceCNs <- readAllOutputBedFiles(dataDir, paramSet, key, numCells, "*__sconce__*.bed") # values from output*bed files
  meanCNs <- readAllOutputBedFiles(dataDir, paramSet, key, numCells, "*__mean.bed")
  eucTrueTree   <- nj(dist(trueCNs, method="euclidean")) # tree built on euclidean distance btn true copy number profiles
  eucSconceTree <- nj(dist(sconceCNs, method="euclidean")) # tree built on euclidean distance btn sconce copy number profiles
  eucMeanTree   <- nj(dist(meanCNs, method="euclidean")) # tree built on euclidean distance btn mean copy number profiles

  #list(eucTrueTree=eucTrueTree, eucSconceTree=eucSconceTree, eucMeanTree=eucMeanTree, eucMedianTree=eucMedianTree, eucModeTree=eucModeTree, cnpTrueTree=cnpTrueTree, cnpSconceTree=cnpSconceTree, cnpMeanTree=cnpMeanTree, cnpMedianTree=cnpMedianTree, cnpModeTree=cnpModeTree, trueDistTree=trueDistTree, inferredDistTree=inferredDistTree)
  list(eucTrueTree=eucTrueTree, eucSconceTree=eucSconceTree, eucMeanTree=eucMeanTree, cnpTrueTree=cnpTrueTree, cnpSconceTree=cnpSconceTree, cnpMeanTree=cnpMeanTree, zzsTrueTree=zzsTrueTree, zzsSconceTree=zzsSconceTree, zzsMeanTree=zzsMeanTree, trueDistTree=trueDistTree, inferredDistTree_t2_t3=inferredDistTree_t2_t3,inferredDistTree_t1=inferredDistTree_t1)
}

calcRFdists <- function(treeList) {
  dists <- data.frame(eucTrue=  treedist(treeList[["trueDistTree"]], treeList[["eucTrueTree"]]),
                      eucSconce=treedist(treeList[["trueDistTree"]], treeList[["eucSconceTree"]]),
                      eucMean=  treedist(treeList[["trueDistTree"]], treeList[["eucMeanTree"]]),
                      eucMedian=treedist(treeList[["trueDistTree"]], treeList[["eucMedianTree"]]),
                      eucMode=  treedist(treeList[["trueDistTree"]], treeList[["eucModeTree"]]),
                      cnpTrue=  treedist(treeList[["trueDistTree"]], treeList[["cnpTrueTree"]]),
                      cnpSconce=treedist(treeList[["trueDistTree"]], treeList[["cnpSconceTree"]]),
                      cnpMean=  treedist(treeList[["trueDistTree"]], treeList[["cnpMeanTree"]]),
                      cnpMedian=treedist(treeList[["trueDistTree"]], treeList[["cnpMedianTree"]]),
                      cnpMode=  treedist(treeList[["trueDistTree"]], treeList[["cnpModeTree"]]),
                      zzsTrue=  treedist(treeList[["trueDistTree"]], treeList[["zzsTrueTree"]]),
                      zzsSconce=treedist(treeList[["trueDistTree"]], treeList[["zzsSconceTree"]]),
                      zzsMean=  treedist(treeList[["trueDistTree"]], treeList[["zzsMeanTree"]]),
                      zzsMedian=treedist(treeList[["trueDistTree"]], treeList[["zzsMedianTree"]]),
                      zzsMode=  treedist(treeList[["trueDistTree"]], treeList[["zzsModeTree"]]),
                      t2_t3=    treedist(treeList[["trueDistTree"]], treeList[["inferredDistTree"]]))
  dists$distMetric <- rownames(dists)
  dists_m <- melt(dists)
  dists_m$cleanTree <- factor(sapply(dists_m$variable, cleanTreeName), levels=orderedTreeNames)

  # reformat for better plotting labels
  dists_m$metric <- sapply(dists_m$cleanTree, FUN=function(x) {s <- unlist(str_split(x, " ")); paste0(s[1], "~", s[2])})
  dists_m$metric <- gsub("dist", "distance", dists_m$metric)
  dists_m$metric <- gsub("inferred~t2+t3", "t[2]+t[3]", dists_m$metric, fixed=T)
  dists_m$metric <- factor(dists_m$metric, levels=metricStrings)

  dists_m$bed <- sapply(dists_m$cleanTree, FUN=function(x) {s <- unlist(str_split(x, " "));s[4]})
  dists_m$bed[is.na(dists_m$bed)] <- "t[2]+t[3]"
  dists_m$bed <- factor(dists_m$bed, levels=bedStrings)

  dists_m
}

calcRFdistsCompareMuts <- function(treeList) {
  dists <- data.frame(eucTrue=  treedist(treeList[["trueDistTree"]], treeList[["eucTrueTree"]]),
                      eucSconce=treedist(treeList[["trueDistTree"]], treeList[["eucSconceTree"]]),
                      eucMean=  treedist(treeList[["trueDistTree"]], treeList[["eucMeanTree"]]),
                      cnpTrue=  treedist(treeList[["trueDistTree"]], treeList[["cnpTrueTree"]]),
                      cnpSconce=treedist(treeList[["trueDistTree"]], treeList[["cnpSconceTree"]]),
                      cnpMean=  treedist(treeList[["trueDistTree"]], treeList[["cnpMeanTree"]]),
                      zzsTrue=  treedist(treeList[["trueDistTree"]], treeList[["zzsTrueTree"]]),
                      zzsSconce=treedist(treeList[["trueDistTree"]], treeList[["zzsSconceTree"]]),
                      zzsMean=  treedist(treeList[["trueDistTree"]], treeList[["zzsMeanTree"]]),
                      t1=       treedist(treeList[["trueDistTree"]], treeList[["inferredDistTree_t1"]]),
                      t2_t3=    treedist(treeList[["trueDistTree"]], treeList[["inferredDistTree_t2_t3"]]))
  dists$distMetric <- rownames(dists)
  dists_m <- melt(dists)
  dists_m$cleanTree <- factor(sapply(dists_m$variable, cleanTreeName), levels=orderedTreeNamesMuts)

  # reformat for better plotting labels
  dists_m$metric <- sapply(dists_m$cleanTree, FUN=function(x) {s <- unlist(str_split(x, " ")); paste0(s[1], "\n", s[2])})
  dists_m$metric <- gsub("dist", "distance", dists_m$metric)
  dists_m$metric <- gsub("inferred\nt2+t3", "t[2]+t[3]", dists_m$metric, fixed=T)
  dists_m$metric <- gsub("inferred\nt1", "dist(t[1])", dists_m$metric, fixed=T)
  dists_m$metric <- factor(dists_m$metric, levels=metricStringsMuts)

  dists_m$bed <- sapply(dists_m$cleanTree, FUN=function(x) {s <- unlist(str_split(x, " "));s[4]})
  dists_m$bed[is.na(dists_m$bed)] <- levels(dists_m$metric)[(dists_m$metric[is.na(dists_m$bed)])]
  dists_m$bed <- factor(gsub("\\]", "", gsub("\\[", "", dists_m$bed)), levels=bedStringsMuts)

  dists_m
}

makeNJtreePlots <- function(treeList, outputFile) {
  njTreesPlotList <- list()
  for(treeName in names(treeList)) {
    plotTitle <- cleanTreeName(treeName) # TODO probably want to reorder and reduce, but not sure how yet
    p <- ggtree(treeList[[treeName]]) + geom_tiplab() + labs(plotTitle=plotTitle)
    njTreesPlotList[[treeName]] <- p
  }
  toSave <- plot_grid(plotlist=njTreesPlotList, align='vh', labels="AUTO")
  png(outputFile, width=9, height=6, units="in", res=600)
  plot(toSave); dev.off()
  njTreesPlotList
}

makeRFdistPlot <- function(dists, plotTitle) {
  metricColors <- hue_pal()(length(metricStrings)) # euc, cnp2cnp, medicc/zzs, t2+t3
  p <- ggplot(subset(dists, distMetric == "symmetric.difference"), aes(x=variable, colour=metric, y=value)) + geom_boxplot() + theme_bw() + theme(legend.title=element_blank(), axis.text.x=element_blank()) + scale_colour_manual(labels=scales::parse_format(), values=metricColors)
  if(!is.na(plotTitle)) {
    p <- p + labs(x="distance metric", y="RF distance", title=plotTitle)
  } else {
    p <- p + labs(x="distance metric", y="RF distance")
  }
  p
}
makeRFdistPlotMuts <- function(dists, plotTitle) {
  metricColors <- hue_pal()(length(metricStringsMuts)) # euc, cnp2cnp, medicc/zzs, t2+t3, t1
  p <- ggplot(subset(dists, distMetric == "symmetric.difference"), aes(x=variable, colour=program, y=value)) + geom_boxplot() + theme_bw() + theme(legend.title=element_blank(), axis.text.x=element_blank()) + scale_colour_manual(labels=scales::parse_format(), values=metricColors)
  if(!is.na(plotTitle)) {
    p <- p + labs(x="distance metric", y="RF distance", title=plotTitle)
  } else {
    p <- p + labs(x="distance metric", y="RF distance")
  }
  p
}
makeRFdistPlotCompareMuts <- function(dists, plotTitle) {
  #p <- ggplot(subset(dists, distMetric == "symmetric.difference"), aes(x=bed, colour=program, y=value)) + geom_boxplot() + theme_bw() + theme(legend.title=element_blank(), axis.title.x=element_blank(), axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) + facet_grid(~metric, space="free", scales="free_x", labeller=label_parsed)
  p <- ggplot(subset(dists, distMetric == "symmetric.difference"), aes(x=bed, fill=program, y=value)) + geom_bar(stat="identity", position="dodge") + theme_bw() + theme(legend.title=element_blank(), axis.title.x=element_blank(), axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) + facet_grid(~metric, space="free", scales="free_x", labeller=label_parsed) + scale_fill_manual(values=sconce2mutColors)
  if(!is.na(plotTitle)) {
    p <- p + labs(x="distance metric", y="RF distance", title=plotTitle)
  } else {
    p <- p + labs(x="distance metric", y="RF distance")
  }
  p
}
# faceted plot: sections/facets are distance metric, true/sconce/mean/median/mode[/t2+t3] are box plots
makeFacetedRFdistPlot <- function(dists, plotTitle) {
  bedColors <- hue_pal()(length(bedStrings)) # sconce, mean, median, mode, true, t2+t3
  p <- ggplot(subset(dists, distMetric == "symmetric.difference"), aes(x=bed, colour=bed, y=value)) + geom_boxplot() + theme_bw() + theme(legend.title=element_blank(), axis.text.x=element_blank()) + facet_grid(~metric, scales="free_x", space="free", labeller=label_parsed) + scale_colour_manual(labels=scales::parse_format(), values=bedColors)

  if(!is.na(plotTitle)) {
    p <- p  + labs(x="distance metric", y="RF distance", title=plotTitle)
  } else {
    p <- p  + labs(x="distance metric", y="RF distance")
  }
  p
}

# write median RF distances to text and tex files for easy copy/paste into latex tables
writeMedianRFdistTexFiles <- function(dists, medianRFdistFile) {
  programs <- levels(dists$variable)
  currParamSets <- unique(dists$paramSet)

  medianDists <- as.data.frame(t(sapply(currParamSets, FUN=function(set) {
    sapply(programs, FUN=function(program) {
      median(subset(dists, paramSet == set & variable == program & distMetric == "symmetric.difference")$value)
    })
  })))
  medianDists$paramSet <- rownames(medianDists)
  medianDists_m <- melt(medianDists)
  medianDists_m$cleanTree <- sapply(medianDists_m$variable, cleanTreeName)

  medianDists_m$metric <- sapply(medianDists_m$cleanTree, FUN=function(x) {s <- unlist(str_split(x, " ")); paste0(s[1], " ", s[2])})
  medianDists_m$metric <- gsub("dist", "distance", medianDists_m$metric)
  medianDists_m$metric <- gsub("inferred ", "", medianDists_m$metric)
  medianDists_m$metric <- factor(medianDists_m$metric, levels=c("Euclidean distance", "cnp2cnp distance", "MEDICC distance", "t2+t3"))

  medianDists_m$bed <- sapply(medianDists_m$cleanTree, FUN=function(x) {s <- unlist(str_split(x, " "));s[4]})
  medianDists_m$bed[is.na(medianDists_m$bed)] <- "t2+t3"
  medianDists_m$bed <- factor(medianDists_m$bed, levels=c("SCONCE", "mean", "median", "mode", "true", "t2+t3"))

  medianDists_w <- dcast(medianDists_m[,c("paramSet", "value", "metric", "bed")],paramSet + metric ~ bed)
  treeNames <- LETTERS[1:length(currParamSets)]
  names(treeNames) <- currParamSets
  medianDists_w$treeName <- treeNames[medianDists_w$paramSet]
  medianDists_w <- medianDists_w[,c("treeName", "metric", "SCONCE", "mean", "median", "mode", "true", "t2+t3", "paramSet")]

  write.table(medianDists_w, file=paste0(medianRFdistFile, "_medianDists.txt"), sep="\t", quote=F, row.names=T, col.names=NA)

  #texColnames <- c("\\textbf{SCONCE}", "\\textbf{one pair}", "\\textbf{mean}", "\\textbf{median}", "\\textbf{mode}", "\\textbf{AneuFinder}") # may need to insert line breaks so the tables fit nicely. see plotJointBreakpointComparison.R for syntax
  texColnames <- paste0("\\textbf{", colnames(medianDists_w), "}")
  names(texColnames) <- colnames(medianDists_w)

  medianDists_w_rounded <- as.data.frame(do.call(cbind, lapply(colnames(medianDists_w), FUN=function(col) {
    if(is.numeric(medianDists_w[,col])) {
      round(medianDists_w[,col], digits=4)
    } else {as.character(medianDists_w[,col])}
  })))
  colnames(medianDists_w_rounded) <- colnames(medianDists_w)

  medianDistsTex <- rbind(texColnames, medianDists_w_rounded)
  medianDistsTex$paramSet <- paste0(" \\\\ \\hline % ", medianDistsTex$paramSet)
  write.table(paste0(apply(medianDistsTex[,-ncol(medianDistsTex)], 1, FUN=function(x) {paste0(x, collapse=" & ")}),  medianDistsTex[,ncol(medianDistsTex)]), file=paste0(medianRFdistFile, "_medianDists.tex"), col.names=F, row.names=F, quote=F)
}


cleanTreeName <- function(f) {
  if(grepl("eucTrue", f)) {
    return("Euclidean dist on true CNPs")
  } else if(grepl("t1", f)) {
    return ("inferred t1")
  } else if(grepl("t2_t3", f) || grepl("inferred", f)) {
    return("inferred t2+t3")
  } else if(grepl("eucSconce", f)) {
    return("Euclidean dist on SCONCE CNPs")
  } else if(grepl("eucMean", f)) {
    return("Euclidean dist on mean CNPs")
  } else if(grepl("eucMedian", f)) {
    return("Euclidean dist on median CNPs")
  } else if(grepl("eucMode", f)) {
    return("Euclidean dist on mode CNPs")
  } else if(grepl("cnpTrue", f)) {
    return("cnp2cnp dist on true CNPs")
  } else if(grepl("cnpSconce", f)) {
    return("cnp2cnp dist on SCONCE CNPs")
  } else if(grepl("cnpMean", f)) {
    return("cnp2cnp dist on mean CNPs")
  } else if(grepl("cnpMedian", f)) {
    return("cnp2cnp dist on median CNPs")
  } else if(grepl("cnpMode", f)) {
    return("cnp2cnp dist on mode CNPs")
  } else if(grepl("zzsTrue", f)) {
    return("MEDICC dist on true CNPs")
  } else if(grepl("zzsSconce", f)) {
    return("MEDICC dist on SCONCE CNPs")
  } else if(grepl("zzsMean", f)) {
    return("MEDICC dist on mean CNPs")
  } else if(grepl("zzsMedian", f)) {
    return("MEDICC dist on median CNPs")
  } else if(grepl("zzsMode", f)) {
    return("MEDICC dist on mode CNPs")
  }
  return(f)
}

orderedTreeNames <- c("Euclidean dist on true CNPs",
                      "Euclidean dist on SCONCE CNPs",
                      "Euclidean dist on mean CNPs",
                      "Euclidean dist on median CNPs",
                      "Euclidean dist on mode CNPs",
                      "cnp2cnp dist on true CNPs",
                      "cnp2cnp dist on SCONCE CNPs",
                      "cnp2cnp dist on mean CNPs",
                      "cnp2cnp dist on median CNPs",
                      "cnp2cnp dist on mode CNPs",
                      "MEDICC dist on true CNPs",
                      "MEDICC dist on SCONCE CNPs",
                      "MEDICC dist on mean CNPs",
                      "MEDICC dist on median CNPs",
                      "MEDICC dist on mode CNPs",
                      "inferred t2+t3")

metricStrings <- c("Euclidean~distance", "cnp2cnp~distance", "MEDICC~distance", "t[2]+t[3]")
bedStrings <- c("SCONCE", "mean", "median", "mode", "true", "t[2]+t[3]")

orderedTreeNamesMuts <- c("Euclidean dist on true CNPs",
                          "Euclidean dist on SCONCE CNPs",
                          "Euclidean dist on mean CNPs",
                          "Euclidean dist on median CNPs",
                          "Euclidean dist on mode CNPs",
                          "cnp2cnp dist on true CNPs",
                          "cnp2cnp dist on SCONCE CNPs",
                          "cnp2cnp dist on mean CNPs",
                          "cnp2cnp dist on median CNPs",
                          "cnp2cnp dist on mode CNPs",
                          "MEDICC dist on true CNPs",
                          "MEDICC dist on SCONCE CNPs",
                          "MEDICC dist on mean CNPs",
                          "MEDICC dist on median CNPs",
                          "MEDICC dist on mode CNPs",
                          "inferred t1",
                          "inferred t2+t3")

metricStringsMuts <- c("Euclidean\ndistance", "cnp2cnp\ndistance", "MEDICC\ndistance", "dist(t[1])", "t[2]+t[3]")
#bedStringsMuts <- c("SCONCE", "mean", "true", "t[1]", "t[2]+t[3]")
bedStringsMuts <- c("SCONCE", "mean", "true", "dist(t1)", "t2+t3")


###########################################
# Tree branch correlation plotting functions
###########################################
# expects results from getTreeBranches
makeCorrPlot <- function(merged_df, plotTitle) {
  # left | right | variable | cell0 | cell1 | inferred | true
  # based on https://stackoverflow.com/questions/7549694/add-regression-line-equation-and-r2-on-graph
  lm_eqn <- function(df){
    m <- lm(inferred ~ true, df)
    rSq <- substitute(italic(r)^2~"="~r2, list(r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(rSq))
  }

  merged_df$variable <- factor(merged_df$variable, levels=c("t1", "t2", "t3", "t2_t3"))
  merged_df$labels <- factor(merged_df$variable)
  levels(merged_df$labels) <- c("t[1]", "t[2]", "t[3]", "t[2]+t[3]")

  lmStr <- ddply(merged_df,.(labels),lm_eqn)

  p <- ggplot(merged_df, aes(x=true, y=inferred, colour=labels)) + geom_point(alpha=0.5) + facet_grid(.~labels, labeller=label_parsed, scales="free") + scale_color_discrete(breaks=levels(merged_df$labels), labels=parse(text=levels(merged_df$labels))) + geom_smooth(method="lm", colour="gray35", show.legend=F) + geom_text(data=lmStr, aes(x=-Inf, y=Inf, label=V1, hjust=0, vjust=1), parse=TRUE, colour="black", show.legend=F) + theme_bw() + theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust=0.5)) + guides(colour=guide_legend(override.aes=list(alpha=1)))
  if(!is.na(plotTitle)) {
    p <- p + labs(colour="branch", x="true", y="inferred", title=plotTitle) 
  } else {
    p <- p + labs(colour="branch", x="true", y="inferred")
  }
  p
}

makeCorrPlotCompareMuts <- function(merged_df, plotTitle) {
  # left | right | variable | cell0 | cell1 | inferred | true
  # based on https://stackoverflow.com/questions/7549694/add-regression-line-equation-and-r2-on-graph
  lm_eqn <- function(df){
    m <- lm(inferred ~ true, df)
    rSq <- substitute(italic(r)^2~"="~r2, list(r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(rSq))
  }

  merged_df$variable <- factor(merged_df$variable, levels=c("t1", "t2", "t3", "t2_t3"))
  merged_df$labels <- factor(merged_df$variable)
  levels(merged_df$labels) <- c("t[1]", "t[2]", "t[3]", "t[2]+t[3]")

  lmStr <- do.call(rbind, lapply(unique(merged_df$program), FUN=function(prog) {
    do.call(rbind, lapply(unique(merged_df$variable), FUN=function(var) {
      data.frame(program=prog, variable=var, text=lm_eqn(subset(merged_df, program == prog & variable == var)))
    }))
  }))
  lmStr$variable <- factor(lmStr$variable, levels=c("t1", "t2", "t3", "t2_t3"))
  lmStr$labels <- factor(lmStr$variable)
  levels(lmStr$labels) <- c("t[1]", "t[2]", "t[3]", "t[2]+t[3]")

  p <- ggplot(merged_df, aes(x=true, y=inferred, colour=labels)) + geom_point(alpha=0.5) + facet_grid(program~labels, labeller=label_parsed, scales="free") + scale_color_discrete(breaks=levels(merged_df$labels), labels=parse(text=levels(merged_df$labels))) + geom_smooth(method="lm", colour="gray35", show.legend=F) + geom_text(data=lmStr, aes(x=-Inf, y=Inf, label=text, hjust=0, vjust=1), parse=TRUE, colour="black", show.legend=F) + theme_bw() + theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust=0.5), legend.position="none") + guides(colour=guide_legend(override.aes=list(alpha=1)))
  if(!is.na(plotTitle)) {
    p <- p + labs(colour="branch", x="true branch length", y="inferred branch length", title=plotTitle) 
  } else {
    p <- p + labs(colour="branch", x="true branch length", y="inferred branch length")
  }

  # shade side strips to match SSE plot coloring
  sourceStripFills <- sconce2mutColors # sconce2 + sconcemut
  sourceTextCols <- c("gray10", "gray10")
  side <- "r"
  g <- ggplot_gtable(ggplot_build(p))
  stripr <- which(grepl(paste0("strip-", side), g$layout$name))
  k <- 1
  for (i in stripr) {
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- sourceStripFills[k]
    j <- which(grepl('title', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$children[[1]]$gp$col <- sourceTextCols[k] # use str(g$grobs) to find elements to change
    k <- k+1
  }
  g
}

###########################################
# newick strings for simulations
###########################################
paramSets <- c("segments/params21", "segments/params22", "segments/params23", "segments/params24", "segments/params25", "segments/params26", "segments/params27", "segments/params28", "segments/params29", "segments/params30", "muts/params14", "muts/params15", "muts/params16", "muts/params17", "muts/params18", "muts/params19", "muts/params20", "muts/params21", "muts/params22", "muts/params23", "muts/params24", "muts/params25", "muts/params26", "muts/params27", "muts/params28", "muts/params29", "muts/params30", "muts/params31", "muts/params32", "muts/params33", "muts/params34")
sp21str <- "(((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((1:0.0078125,2:0.0078125):0.0078125,3:0.015625):0.0078125,4:0.0234375):0.0078125,5:0.03125):0.0078125,6:0.0390625):0.0078125,7:0.046875):0.0078125,8:0.0546875):0.0078125,9:0.0625):0.0078125,10:0.0703125):0.0078125,11:0.078125):0.0078125,12:0.0859375):0.0078125,13:0.09375):0.0078125,14:0.1015625):0.0078125,15:0.109375):0.0078125,16:0.1171875):0.0078125,17:0.125):0.0078125,18:0.1328125):0.0078125,19:0.140625):0.0078125,20:0.1484375):0.0078125,21:0.15625):0.0078125,22:0.1640625):0.0078125,23:0.171875):0.0078125,24:0.1796875):0.0078125,25:0.1875):0.0078125,26:0.1953125):0.0078125,27:0.203125):0.0078125,28:0.2109375):0.0078125,29:0.21875):0.0078125,30:0.2265625):0.0078125,31:0.234375):0.0078125,32:0.2421875):0.0078125,33:0.25):0.0078125,34:0.2578125):0.0078125,35:0.265625):0.0078125,36:0.2734375):0.0078125,37:0.28125):0.0078125,38:0.2890625):0.0078125,39:0.296875):0.0078125,40:0.3046875):0.0078125,41:0.3125):0.0078125,42:0.3203125):0.0078125,43:0.328125):0.0078125,44:0.3359375):0.0078125,45:0.34375):0.0078125,46:0.3515625):0.0078125,47:0.359375):0.0078125,48:0.3671875):0.0078125,49:0.375):0.0078125,50:0.3828125):0.0078125,51:0.390625):0.0078125,52:0.3984375):0.0078125,53:0.40625):0.0078125,54:0.4140625):0.0078125,55:0.421875):0.0078125,56:0.4296875):0.0078125,57:0.4375):0.0078125,58:0.4453125):0.0078125,59:0.453125):0.0078125,60:0.4609375):0.0078125,61:0.46875):0.0078125,62:0.4765625):0.0078125,63:0.484375):0.0078125,64:0.4921875):0.0078125,65:0.5):0.0078125,66:0.5078125):0.0078125,67:0.515625):0.0078125,68:0.5234375):0.0078125,69:0.53125):0.0078125,70:0.5390625):0.0078125,71:0.546875):0.0078125,72:0.5546875):0.0078125,73:0.5625):0.0078125,74:0.5703125):0.0078125,75:0.578125):0.0078125,76:0.5859375):0.0078125,77:0.59375):0.0078125,78:0.6015625):0.0078125,79:0.609375):0.0078125,80:0.6171875):0.0078125,81:0.625):0.0078125,82:0.6328125):0.0078125,83:0.640625):0.0078125,84:0.6484375):0.0078125,85:0.65625):0.0078125,86:0.6640625):0.0078125,87:0.671875):0.0078125,88:0.6796875):0.0078125,89:0.6875):0.0078125,90:0.6953125):0.0078125,91:0.703125):0.0078125,92:0.7109375):0.0078125,93:0.71875):0.0078125,94:0.7265625):0.0078125,95:0.734375):0.0078125,96:0.7421875):0.0078125,97:0.75):0.0078125,98:0.7578125):0.0078125,99:0.765625):0.0078125,100:0.7734375):0.0078125,101:0.78125):0.0078125,102:0.7890625):0.0078125,103:0.796875):0.0078125,104:0.8046875):0.0078125,105:0.8125):0.0078125,106:0.8203125):0.0078125,107:0.828125):0.0078125,108:0.8359375):0.0078125,109:0.84375):0.0078125,110:0.8515625):0.0078125,111:0.859375):0.0078125,112:0.8671875):0.0078125,113:0.875):0.0078125,114:0.8828125):0.0078125,115:0.890625):0.0078125,116:0.8984375):0.0078125,117:0.90625):0.0078125,118:0.9140625):0.0078125,119:0.921875):0.0078125,120:0.9296875):0.0078125,121:0.9375):0.0078125,122:0.9453125):0.0078125,123:0.953125):0.0078125,124:0.9609375):0.0078125,125:0.96875):0.0078125,126:0.9765625):0.0078125,127:0.984375):0.0078125,128:0.9921875);"
sp22str <- sp27str <- sp28str <- "(((((((1:0.125,2:0.125):0.125,(3:0.125,4:0.125):0.125):0.125,((5:0.125,6:0.125):0.125,(7:0.125,8:0.125):0.125):0.125):0.125,(((9:0.125,10:0.125):0.125,(11:0.125,12:0.125):0.125):0.125,((13:0.125,14:0.125):0.125,(15:0.125,16:0.125):0.125):0.125):0.125):0.125,((((17:0.125,18:0.125):0.125,(19:0.125,20:0.125):0.125):0.125,((21:0.125,22:0.125):0.125,(23:0.125,24:0.125):0.125):0.125):0.125,(((25:0.125,26:0.125):0.125,(27:0.125,28:0.125):0.125):0.125,((29:0.125,30:0.125):0.125,(31:0.125,32:0.125):0.125):0.125):0.125):0.125):0.125,(((((33:0.125,34:0.125):0.125,(35:0.125,36:0.125):0.125):0.125,((37:0.125,38:0.125):0.125,(39:0.125,40:0.125):0.125):0.125):0.125,(((41:0.125,42:0.125):0.125,(43:0.125,44:0.125):0.125):0.125,((45:0.125,46:0.125):0.125,(47:0.125,48:0.125):0.125):0.125):0.125):0.125,((((49:0.125,50:0.125):0.125,(51:0.125,52:0.125):0.125):0.125,((53:0.125,54:0.125):0.125,(55:0.125,56:0.125):0.125):0.125):0.125,(((57:0.125,58:0.125):0.125,(59:0.125,60:0.125):0.125):0.125,((61:0.125,62:0.125):0.125,(63:0.125,64:0.125):0.125):0.125):0.125):0.125):0.125):0.125,((((((65:0.125,66:0.125):0.125,(67:0.125,68:0.125):0.125):0.125,((69:0.125,70:0.125):0.125,(71:0.125,72:0.125):0.125):0.125):0.125,(((73:0.125,74:0.125):0.125,(75:0.125,76:0.125):0.125):0.125,((77:0.125,78:0.125):0.125,(79:0.125,80:0.125):0.125):0.125):0.125):0.125,((((81:0.125,82:0.125):0.125,(83:0.125,84:0.125):0.125):0.125,((85:0.125,86:0.125):0.125,(87:0.125,88:0.125):0.125):0.125):0.125,(((89:0.125,90:0.125):0.125,(91:0.125,92:0.125):0.125):0.125,((93:0.125,94:0.125):0.125,(95:0.125,96:0.125):0.125):0.125):0.125):0.125):0.125,(((((97:0.125,98:0.125):0.125,(99:0.125,100:0.125):0.125):0.125,((101:0.125,102:0.125):0.125,(103:0.125,104:0.125):0.125):0.125):0.125,(((105:0.125,106:0.125):0.125,(107:0.125,108:0.125):0.125):0.125,((109:0.125,110:0.125):0.125,(111:0.125,112:0.125):0.125):0.125):0.125):0.125,((((113:0.125,114:0.125):0.125,(115:0.125,116:0.125):0.125):0.125,((117:0.125,118:0.125):0.125,(119:0.125,120:0.125):0.125):0.125):0.125,(((121:0.125,122:0.125):0.125,(123:0.125,124:0.125):0.125):0.125,((125:0.125,126:0.125):0.125,(127:0.125,128:0.125):0.125):0.125):0.125):0.125):0.125):0.125);"
sp23str <- sp25str <- "(((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((1:0.0078125,2:0.0078125):0.0078125,3:0.0078125):0.0078125,4:0.0078125):0.0078125,5:0.0078125):0.0078125,6:0.0078125):0.0078125,7:0.0078125):0.0078125,8:0.0078125):0.0078125,9:0.0078125):0.0078125,10:0.0078125):0.0078125,11:0.0078125):0.0078125,12:0.0078125):0.0078125,13:0.0078125):0.0078125,14:0.0078125):0.0078125,15:0.0078125):0.0078125,16:0.0078125):0.0078125,17:0.0078125):0.0078125,18:0.0078125):0.0078125,19:0.0078125):0.0078125,20:0.0078125):0.0078125,21:0.0078125):0.0078125,22:0.0078125):0.0078125,23:0.0078125):0.0078125,24:0.0078125):0.0078125,25:0.0078125):0.0078125,26:0.0078125):0.0078125,27:0.0078125):0.0078125,28:0.0078125):0.0078125,29:0.0078125):0.0078125,30:0.0078125):0.0078125,31:0.0078125):0.0078125,32:0.0078125):0.0078125,33:0.0078125):0.0078125,34:0.0078125):0.0078125,35:0.0078125):0.0078125,36:0.0078125):0.0078125,37:0.0078125):0.0078125,38:0.0078125):0.0078125,39:0.0078125):0.0078125,40:0.0078125):0.0078125,41:0.0078125):0.0078125,42:0.0078125):0.0078125,43:0.0078125):0.0078125,44:0.0078125):0.0078125,45:0.0078125):0.0078125,46:0.0078125):0.0078125,47:0.0078125):0.0078125,48:0.0078125):0.0078125,49:0.0078125):0.0078125,50:0.0078125):0.0078125,51:0.0078125):0.0078125,52:0.0078125):0.0078125,53:0.0078125):0.0078125,54:0.0078125):0.0078125,55:0.0078125):0.0078125,56:0.0078125):0.0078125,57:0.0078125):0.0078125,58:0.0078125):0.0078125,59:0.0078125):0.0078125,60:0.0078125):0.0078125,61:0.0078125):0.0078125,62:0.0078125):0.0078125,63:0.0078125):0.0078125,64:0.0078125):0.0078125,65:0.0078125):0.0078125,66:0.0078125):0.0078125,67:0.0078125):0.0078125,68:0.0078125):0.0078125,69:0.0078125):0.0078125,70:0.0078125):0.0078125,71:0.0078125):0.0078125,72:0.0078125):0.0078125,73:0.0078125):0.0078125,74:0.0078125):0.0078125,75:0.0078125):0.0078125,76:0.0078125):0.0078125,77:0.0078125):0.0078125,78:0.0078125):0.0078125,79:0.0078125):0.0078125,80:0.0078125):0.0078125,81:0.0078125):0.0078125,82:0.0078125):0.0078125,83:0.0078125):0.0078125,84:0.0078125):0.0078125,85:0.0078125):0.0078125,86:0.0078125):0.0078125,87:0.0078125):0.0078125,88:0.0078125):0.0078125,89:0.0078125):0.0078125,90:0.0078125):0.0078125,91:0.0078125):0.0078125,92:0.0078125):0.0078125,93:0.0078125):0.0078125,94:0.0078125):0.0078125,95:0.0078125):0.0078125,96:0.0078125):0.0078125,97:0.0078125):0.0078125,98:0.0078125):0.0078125,99:0.0078125):0.0078125,100:0.0078125):0.0078125,101:0.0078125):0.0078125,102:0.0078125):0.0078125,103:0.0078125):0.0078125,104:0.0078125):0.0078125,105:0.0078125):0.0078125,106:0.0078125):0.0078125,107:0.0078125):0.0078125,108:0.0078125):0.0078125,109:0.0078125):0.0078125,110:0.0078125):0.0078125,111:0.0078125):0.0078125,112:0.0078125):0.0078125,113:0.0078125):0.0078125,114:0.0078125):0.0078125,115:0.0078125):0.0078125,116:0.0078125):0.0078125,117:0.0078125):0.0078125,118:0.0078125):0.0078125,119:0.0078125):0.0078125,120:0.0078125):0.0078125,121:0.0078125):0.0078125,122:0.0078125):0.0078125,123:0.0078125):0.0078125,124:0.0078125):0.0078125,125:0.0078125):0.0078125,126:0.0078125):0.0078125,127:0.0078125):0.0078125,128:0.0078125);"
sp24str <- "(((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((1:0.015625,2:0.015625):0.0078125,3:0.015625):0.0078125,4:0.015625):0.0078125,5:0.015625):0.0078125,6:0.015625):0.0078125,7:0.015625):0.0078125,8:0.015625):0.0078125,9:0.015625):0.0078125,10:0.015625):0.0078125,11:0.015625):0.0078125,12:0.015625):0.0078125,13:0.015625):0.0078125,14:0.015625):0.0078125,15:0.015625):0.0078125,16:0.015625):0.0078125,17:0.015625):0.0078125,18:0.015625):0.0078125,19:0.015625):0.0078125,20:0.015625):0.0078125,21:0.015625):0.0078125,22:0.015625):0.0078125,23:0.015625):0.0078125,24:0.015625):0.0078125,25:0.015625):0.0078125,26:0.015625):0.0078125,27:0.015625):0.0078125,28:0.015625):0.0078125,29:0.015625):0.0078125,30:0.015625):0.0078125,31:0.015625):0.0078125,32:0.015625):0.0078125,33:0.015625):0.0078125,34:0.015625):0.0078125,35:0.015625):0.0078125,36:0.015625):0.0078125,37:0.015625):0.0078125,38:0.015625):0.0078125,39:0.015625):0.0078125,40:0.015625):0.0078125,41:0.015625):0.0078125,42:0.015625):0.0078125,43:0.015625):0.0078125,44:0.015625):0.0078125,45:0.015625):0.0078125,46:0.015625):0.0078125,47:0.015625):0.0078125,48:0.015625):0.0078125,49:0.015625):0.0078125,50:0.015625):0.0078125,51:0.015625):0.0078125,52:0.015625):0.0078125,53:0.015625):0.0078125,54:0.015625):0.0078125,55:0.015625):0.0078125,56:0.015625):0.0078125,57:0.015625):0.0078125,58:0.015625):0.0078125,59:0.015625):0.0078125,60:0.015625):0.0078125,61:0.015625):0.0078125,62:0.015625):0.0078125,63:0.015625):0.0078125,64:0.015625):0.0078125,65:0.015625):0.0078125,66:0.015625):0.0078125,67:0.015625):0.0078125,68:0.015625):0.0078125,69:0.015625):0.0078125,70:0.015625):0.0078125,71:0.015625):0.0078125,72:0.015625):0.0078125,73:0.015625):0.0078125,74:0.015625):0.0078125,75:0.015625):0.0078125,76:0.015625):0.0078125,77:0.015625):0.0078125,78:0.015625):0.0078125,79:0.015625):0.0078125,80:0.015625):0.0078125,81:0.015625):0.0078125,82:0.015625):0.0078125,83:0.015625):0.0078125,84:0.015625):0.0078125,85:0.015625):0.0078125,86:0.015625):0.0078125,87:0.015625):0.0078125,88:0.015625):0.0078125,89:0.015625):0.0078125,90:0.015625):0.0078125,91:0.015625):0.0078125,92:0.015625):0.0078125,93:0.015625):0.0078125,94:0.015625):0.0078125,95:0.015625):0.0078125,96:0.015625):0.0078125,97:0.015625):0.0078125,98:0.015625):0.0078125,99:0.015625):0.0078125,100:0.015625):0.0078125,101:0.015625):0.0078125,102:0.015625):0.0078125,103:0.015625):0.0078125,104:0.015625):0.0078125,105:0.015625):0.0078125,106:0.015625):0.0078125,107:0.015625):0.0078125,108:0.015625):0.0078125,109:0.015625):0.0078125,110:0.015625):0.0078125,111:0.015625):0.0078125,112:0.015625):0.0078125,113:0.015625):0.0078125,114:0.015625):0.0078125,115:0.015625):0.0078125,116:0.015625):0.0078125,117:0.015625):0.0078125,118:0.015625):0.0078125,119:0.015625):0.0078125,120:0.015625):0.0078125,121:0.015625):0.0078125,122:0.015625):0.0078125,123:0.015625):0.0078125,124:0.015625):0.0078125,125:0.015625):0.0078125,126:0.015625):0.0078125,127:0.015625):0.0078125,128:0.015625);"
sp26str <- "(((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((1:0.0078125,2:0.0078125):0.0078125,3:0.0116008332014688):0.0078125,4:0.0155916583871406):0.0078125,5:0.0187593545442356):0.0078125,6:0.0213880438928388):0.0078125,7:0.0236355600838684):0.0078125,8:0.0255989540548129):0.0078125,9:0.0273422630494552):0.0078125,10:0.0289100028577532):0.0078125,11:0.030334390020206):0.0078125,12:0.0316395000566371):0.0078125,13:0.0328438022427875):0.0078125,14:0.0339617775837302):0.0078125,15:0.0350049915886093):0.0078125,16:0.0359828286429374):0.0078125,17:0.036903008609003):0.0078125,18:0.0377719588833015):0.0078125,19:0.0385950879118784):0.0078125,20:0.0393769899318549):0.0078125,21:0.0401216007115761):0.0078125,22:0.0408323177288084):0.0078125,23:0.0415120941115399):0.0078125,24:0.0421635129313762):0.0078125,25:0.0427888465850756):0.0078125,26:0.0433901047189515):0.0078125,27:0.0439690732514026):0.0078125,28:0.0445273464075131):0.0078125,29:0.0450663532160049):0.0078125,30:0.0455873795792758):0.0078125,31:0.0460915867756631):0.0078125,32:0.0465800270645798):0.0078125,33:0.0470536569225207):0.0078125,34:0.0475133483289498):0.0078125,35:0.0479598984370744):0.0078125,36:0.0483940378992214):0.0078125,37:0.0488164380653954):0.0078125,38:0.0492277172332503):0.0078125,39:0.049628446095653):0.0078125,40:0.0500191525063924):0.0078125,41:0.0504003256639661):0.0078125,42:0.0507724197966948):0.0078125,43:0.0511358574188424):0.0078125,44:0.051491032216317):0.0078125,45:0.0518383116114016):0.0078125,46:0.0521780390484309):0.0078125,47:0.0525105360360778):0.0078125,48:0.0528361039767072):0.0078125,49:0.0531550258088966):0.0078125,50:0.0534675674855661):0.0078125,51:0.0537739793070762):0.0078125,52:0.0540744971260433):0.0078125,53:0.0543693434384093):0.0078125,54:0.054658728373412):0.0078125,55:0.054942850593497):0.0078125,56:0.0552218981138258):0.0078125,57:0.0554960490498514):0.0078125,58:0.0557654723004074):0.0078125,59:0.0560303281728737):0.0078125,60:0.0562907689562142):0.0078125,61:0.0565469394470188):0.0078125,62:0.0567989774330987):0.0078125,63:0.0570470141386813):0.0078125,64:0.057291174634807):0.0078125,65:0.057531578218141):0.0078125,66:0.0577683387610733):0.0078125,67:0.0580015650356803):0.0078125,68:0.0582313610138545):0.0078125,69:0.0584578261456767):0.0078125,70:0.0586810556178959):0.0078125,71:0.0589011405941983):0.0078125,72:0.0591181684387847):0.0078125,73:0.0593322229246294):0.0078125,74:0.0595433844276601):0.0078125,75:0.0597517301079894):0.0078125,76:0.0599573340792189):0.0078125,77:0.0601602675667444):0.0078125,78:0.0603605990559111):0.0078125,79:0.0605583944307883):0.0078125,80:0.0607537171042681):0.0078125,81:0.0609466281401304):0.0078125,82:0.0611371863676629):0.0078125,83:0.0613254484893743):0.0078125,84:0.061511469182294):0.0078125,85:0.0616953011933111):0.0078125,86:0.0618769954289694):0.0078125,87:0.0620566010400999):0.0078125,88:0.0622341655016433):0.0078125,89:0.0624097346879875):0.0078125,90:0.0625833529441181):0.0078125,91:0.0627550631528585):0.0078125,92:0.0629249067984554):0.0078125,93:0.0630929240267443):0.0078125,94:0.0632591537021137):0.0078125,95:0.0634236334614717):0.0078125,96:0.0635863997654003):0.0078125,97:0.0637474879466715):0.0078125,98:0.0639069322562904):0.0078125,99:0.0640647659072084):0.0078125,100:0.0642210211158543):0.0078125,101:0.0643757291416057):0.0078125,102:0.0645289203243274):0.0078125,103:0.064680624120086):0.0078125,104:0.0648308691351476):0.0078125,105:0.0649796831583561):0.0078125,106:0.0651270931919842):0.0078125,107:0.0652731254811416):0.0078125,108:0.0654178055418231):0.0078125,109:0.0655611581876684):0.0078125,110:0.0657032075555058):0.0078125,111:0.0658439771297455):0.0078125,112:0.0659834897656831):0.0078125,113:0.0661217677117721):0.0078125,114:0.0662588326309198):0.0078125,115:0.0663947056208564):0.0078125,116:0.066529407233626):0.0078125,117:0.0666629574942454):0.0078125,118:0.0667953759185718):0.0078125,119:0.0669266815304186):0.0078125,120:0.0670568928779594):0.0078125,121:0.0671860280494532):0.0078125,122:0.0673141046883257):0.0078125,123:0.0674411400076365):0.0078125,124:0.0675671508039636):0.0078125,125:0.0676921534707327):0.0078125,126:0.0678161640110162):0.0078125,127:0.0679391980498297):0.0078125,128:0.0680612708459476);"
#sp29str <- "(((1:0.05,2:0.05):0.05,(3:0.05,4:0.05):0.05):0.9,((5:0.05,6:0.05):0.05,(7:0.05,8:0.05):0.05):0.9);"
#sp30str <- "(((1:0.9,2:0.9):0.05,(3:0.9,4:0.9):0.05):0.05,((5:0.9,6:0.9):0.05,(7:0.9,8:0.9):0.05):0.05);"
sp29str <- "((((1:0.033,2:0.033):0.033,(3:0.033,4:0.033):0.033):0.033,((5:0.033,6:0.033):0.033,(7:0.033,8:0.033):0.033):0.033):0.9,(((9:0.033,10:0.033):0.033,(11:0.033,12:0.033):0.033):0.033,((13:0.033,14:0.033):0.033,(15:0.033,16:0.033):0.033):0.033):0.9);"
sp30str <- "((((1:0.9,2:0.9):0.033,(3:0.9,4:0.9):0.033):0.033,((5:0.9,6:0.9):0.033,(7:0.9,8:0.9):0.033):0.033):0.033,(((9:0.9,10:0.9):0.033,(11:0.9,12:0.9):0.033):0.033,((13:0.9,14:0.9):0.033,(15:0.9,16:0.9):0.033):0.033):0.033);"
mp14str <- sp21str
mp15str <- sp22str
mp16str <- sp25str
mp17str <- sp26str
mp18str <- mp14str
mp19str <- mp14str
mp20str <- mp14str
mp21str <- mp18str
mp22str <- "(((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((1:0.0078740157480315,2:0.0078740157480315):0.0078740157480315,3:0.015748031496063):0.0078740157480315,4:0.0236220472440945):0.0078740157480315,5:0.031496062992126):0.0078740157480315,6:0.0393700787401575):0.0078740157480315,7:0.047244094488189):0.0078740157480315,8:0.0551181102362205):0.0078740157480315,9:0.062992125984252):0.0078740157480315,10:0.0708661417322835):0.0078740157480315,11:0.078740157480315):0.0078740157480315,12:0.0866141732283465):0.0078740157480315,13:0.094488188976378):0.0078740157480315,14:0.102362204724409):0.0078740157480315,15:0.110236220472441):0.0078740157480315,16:0.118110236220472):0.0078740157480315,17:0.125984251968504):0.0078740157480315,18:0.133858267716535):0.0078740157480315,19:0.141732283464567):0.0078740157480315,20:0.149606299212598):0.0078740157480315,21:0.15748031496063):0.0078740157480315,22:0.165354330708661):0.0078740157480315,23:0.173228346456693):0.0078740157480315,24:0.181102362204724):0.0078740157480315,25:0.188976377952756):0.0078740157480315,26:0.196850393700787):0.0078740157480315,27:0.204724409448819):0.0078740157480315,28:0.21259842519685):0.0078740157480315,29:0.220472440944882):0.0078740157480315,30:0.228346456692913):0.0078740157480315,31:0.236220472440945):0.0078740157480315,32:0.244094488188976):0.0078740157480315,33:0.251968503937008):0.0078740157480315,34:0.259842519685039):0.0078740157480315,35:0.267716535433071):0.0078740157480315,36:0.275590551181102):0.0078740157480315,37:0.283464566929134):0.0078740157480315,38:0.291338582677165):0.0078740157480315,39:0.299212598425197):0.0078740157480315,40:0.307086614173228):0.0078740157480315,41:0.31496062992126):0.0078740157480315,42:0.322834645669291):0.0078740157480315,43:0.330708661417323):0.0078740157480315,44:0.338582677165354):0.0078740157480315,45:0.346456692913386):0.0078740157480315,46:0.354330708661417):0.0078740157480315,47:0.362204724409449):0.0078740157480315,48:0.37007874015748):0.0078740157480315,49:0.377952755905512):0.0078740157480315,50:0.385826771653543):0.0078740157480315,51:0.393700787401575):0.0078740157480315,52:0.401574803149606):0.0078740157480315,53:0.409448818897638):0.0078740157480315,54:0.417322834645669):0.0078740157480315,55:0.425196850393701):0.0078740157480315,56:0.433070866141732):0.0078740157480315,57:0.440944881889764):0.0078740157480315,58:0.448818897637795):0.0078740157480315,59:0.456692913385827):0.0078740157480315,60:0.464566929133858):0.0078740157480315,61:0.47244094488189):0.0078740157480315,62:0.480314960629921):0.0078740157480315,63:0.488188976377953):0.0078740157480315,64:0.496062992125984):0.0078740157480315,65:0.503937007874016):0.0078740157480315,66:0.511811023622047):0.0078740157480315,67:0.519685039370079):0.0078740157480315,68:0.52755905511811):0.0078740157480315,69:0.535433070866142):0.0078740157480315,70:0.543307086614173):0.0078740157480315,71:0.551181102362205):0.0078740157480315,72:0.559055118110236):0.0078740157480315,73:0.566929133858268):0.0078740157480315,74:0.574803149606299):0.0078740157480315,75:0.582677165354331):0.0078740157480315,76:0.590551181102362):0.0078740157480315,77:0.598425196850394):0.0078740157480315,78:0.606299212598425):0.0078740157480315,79:0.614173228346457):0.0078740157480315,80:0.622047244094488):0.0078740157480315,81:0.62992125984252):0.0078740157480315,82:0.637795275590551):0.0078740157480315,83:0.645669291338583):0.0078740157480315,84:0.653543307086614):0.0078740157480315,85:0.661417322834646):0.0078740157480315,86:0.669291338582677):0.0078740157480315,87:0.677165354330709):0.0078740157480315,88:0.68503937007874):0.0078740157480315,89:0.692913385826772):0.0078740157480315,90:0.700787401574803):0.0078740157480315,91:0.708661417322835):0.0078740157480315,92:0.716535433070866):0.0078740157480315,93:0.724409448818898):0.0078740157480315,94:0.732283464566929):0.0078740157480315,95:0.740157480314961):0.0078740157480315,96:0.748031496062992):0.0078740157480315,97:0.755905511811024):0.0078740157480315,98:0.763779527559055):0.0078740157480315,99:0.771653543307087):0.0078740157480315,100:0.779527559055118):0.0078740157480315,101:0.78740157480315):0.0078740157480315,102:0.795275590551181):0.0078740157480315,103:0.803149606299213):0.0078740157480315,104:0.811023622047244):0.0078740157480315,105:0.818897637795276):0.0078740157480315,106:0.826771653543307):0.0078740157480315,107:0.834645669291339):0.0078740157480315,108:0.84251968503937):0.0078740157480315,109:0.850393700787402):0.0078740157480315,110:0.858267716535433):0.0078740157480315,111:0.866141732283465):0.0078740157480315,112:0.874015748031496):0.0078740157480315,113:0.881889763779528):0.0078740157480315,114:0.889763779527559):0.0078740157480315,115:0.897637795275591):0.0078740157480315,116:0.905511811023622):0.0078740157480315,117:0.913385826771653):0.0078740157480315,118:0.921259842519685):0.0078740157480315,119:0.929133858267717):0.0078740157480315,120:0.937007874015748):0.0078740157480315,121:0.94488188976378):0.0078740157480315,122:0.952755905511811):0.0078740157480315,123:0.960629921259842):0.0078740157480315,124:0.968503937007874):0.0078740157480315,125:0.976377952755906):0.0078740157480315,126:0.984251968503937):0.0078740157480315,127:0.992125984251969):0.0078740157480315,128:1);"
mp23str <- "((((1:.25,2:.25):.25,3:.5):.25,4:.75):.25,5:1);"
mp24str <- mp28str <- mp29str <- mp30str <- mp31str <- "(((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((1:0.0078740157480315,2:0.0078740157480315):0.0078740157480315,3:0.015748031496063):0.0078740157480315,4:0.0236220472440945):0.0078740157480315,5:0.031496062992126):0.0078740157480315,6:0.0393700787401575):0.0078740157480315,7:0.047244094488189):0.0078740157480315,8:0.0551181102362205):0.0078740157480315,9:0.062992125984252):0.0078740157480315,10:0.0708661417322835):0.0078740157480315,11:0.078740157480315):0.0078740157480315,12:0.0866141732283465):0.0078740157480315,13:0.094488188976378):0.0078740157480315,14:0.102362204724409):0.0078740157480315,15:0.110236220472441):0.0078740157480315,16:0.118110236220472):0.0078740157480315,17:0.125984251968504):0.0078740157480315,18:0.133858267716535):0.0078740157480315,19:0.141732283464567):0.0078740157480315,20:0.149606299212598):0.0078740157480315,21:0.15748031496063):0.0078740157480315,22:0.165354330708661):0.0078740157480315,23:0.173228346456693):0.0078740157480315,24:0.181102362204724):0.0078740157480315,25:0.188976377952756):0.0078740157480315,26:0.196850393700787):0.0078740157480315,27:0.204724409448819):0.0078740157480315,28:0.21259842519685):0.0078740157480315,29:0.220472440944882):0.0078740157480315,30:0.228346456692913):0.0078740157480315,31:0.236220472440945):0.0078740157480315,32:0.244094488188976):0.0078740157480315,33:0.251968503937008):0.0078740157480315,34:0.259842519685039):0.0078740157480315,35:0.267716535433071):0.0078740157480315,36:0.275590551181102):0.0078740157480315,37:0.283464566929134):0.0078740157480315,38:0.291338582677165):0.0078740157480315,39:0.299212598425197):0.0078740157480315,40:0.307086614173228):0.0078740157480315,41:0.31496062992126):0.0078740157480315,42:0.322834645669291):0.0078740157480315,43:0.330708661417323):0.0078740157480315,44:0.338582677165354):0.0078740157480315,45:0.346456692913386):0.0078740157480315,46:0.354330708661417):0.0078740157480315,47:0.362204724409449):0.0078740157480315,48:0.37007874015748):0.0078740157480315,49:0.377952755905512):0.0078740157480315,50:0.385826771653543):0.0078740157480315,51:0.393700787401575):0.0078740157480315,52:0.401574803149606):0.0078740157480315,53:0.409448818897638):0.0078740157480315,54:0.417322834645669):0.0078740157480315,55:0.425196850393701):0.0078740157480315,56:0.433070866141732):0.0078740157480315,57:0.440944881889764):0.0078740157480315,58:0.448818897637795):0.0078740157480315,59:0.456692913385827):0.0078740157480315,60:0.464566929133858):0.0078740157480315,61:0.47244094488189):0.0078740157480315,62:0.480314960629921):0.0078740157480315,63:0.488188976377953):0.0078740157480315,64:0.496062992125984):0.0078740157480315,65:0.503937007874016):0.0078740157480315,66:0.511811023622047):0.0078740157480315,67:0.519685039370079):0.0078740157480315,68:0.52755905511811):0.0078740157480315,69:0.535433070866142):0.0078740157480315,70:0.543307086614173):0.0078740157480315,71:0.551181102362205):0.0078740157480315,72:0.559055118110236):0.0078740157480315,73:0.566929133858268):0.0078740157480315,74:0.574803149606299):0.0078740157480315,75:0.582677165354331):0.0078740157480315,76:0.590551181102362):0.0078740157480315,77:0.598425196850394):0.0078740157480315,78:0.606299212598425):0.0078740157480315,79:0.614173228346457):0.0078740157480315,80:0.622047244094488):0.0078740157480315,81:0.62992125984252):0.0078740157480315,82:0.637795275590551):0.0078740157480315,83:0.645669291338583):0.0078740157480315,84:0.653543307086614):0.0078740157480315,85:0.661417322834646):0.0078740157480315,86:0.669291338582677):0.0078740157480315,87:0.677165354330709):0.0078740157480315,88:0.68503937007874):0.0078740157480315,89:0.692913385826772):0.0078740157480315,90:0.700787401574803):0.0078740157480315,91:0.708661417322835):0.0078740157480315,92:0.716535433070866):0.0078740157480315,93:0.724409448818898):0.0078740157480315,94:0.732283464566929):0.0078740157480315,95:0.740157480314961):0.0078740157480315,96:0.748031496062992):0.0078740157480315,97:0.755905511811024):0.0078740157480315,98:0.763779527559055):0.0078740157480315,99:0.771653543307087):0.0078740157480315,100:0.779527559055118):0.0078740157480315,101:0.78740157480315):0.0078740157480315,102:0.795275590551181):0.0078740157480315,103:0.803149606299213):0.0078740157480315,104:0.811023622047244):0.0078740157480315,105:0.818897637795276):0.0078740157480315,106:0.826771653543307):0.0078740157480315,107:0.834645669291339):0.0078740157480315,108:0.84251968503937):0.0078740157480315,109:0.850393700787402):0.0078740157480315,110:0.858267716535433):0.0078740157480315,111:0.866141732283465):0.0078740157480315,112:0.874015748031496):0.0078740157480315,113:0.881889763779528):0.0078740157480315,114:0.889763779527559):0.0078740157480315,115:0.897637795275591):0.0078740157480315,116:0.905511811023622):0.0078740157480315,117:0.913385826771653):0.0078740157480315,118:0.921259842519685):0.0078740157480315,119:0.929133858267717):0.0078740157480315,120:0.937007874015748):0.0078740157480315,121:0.94488188976378):0.0078740157480315,122:0.952755905511811):0.0078740157480315,123:0.960629921259842):0.0078740157480315,124:0.968503937007874):0.0078740157480315,125:0.976377952755906):0.0078740157480315,126:0.984251968503937):0.0078740157480315,127:0.992125984251969):0.0078740157480315,128:1);"
mp25str <- mp15str
mp26str <- "(((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((1:0.0078740157480315,2:0.0078740157480315):0.0078740157480315,3:0.0078740157480315):0.0078740157480315,4:0.0078740157480315):0.0078740157480315,5:0.0078740157480315):0.0078740157480315,6:0.0078740157480315):0.0078740157480315,7:0.0078740157480315):0.0078740157480315,8:0.0078740157480315):0.0078740157480315,9:0.0078740157480315):0.0078740157480315,10:0.0078740157480315):0.0078740157480315,11:0.0078740157480315):0.0078740157480315,12:0.0078740157480315):0.0078740157480315,13:0.0078740157480315):0.0078740157480315,14:0.0078740157480315):0.0078740157480315,15:0.0078740157480315):0.0078740157480315,16:0.0078740157480315):0.0078740157480315,17:0.0078740157480315):0.0078740157480315,18:0.0078740157480315):0.0078740157480315,19:0.0078740157480315):0.0078740157480315,20:0.0078740157480315):0.0078740157480315,21:0.0078740157480315):0.0078740157480315,22:0.0078740157480315):0.0078740157480315,23:0.0078740157480315):0.0078740157480315,24:0.0078740157480315):0.0078740157480315,25:0.0078740157480315):0.0078740157480315,26:0.0078740157480315):0.0078740157480315,27:0.0078740157480315):0.0078740157480315,28:0.0078740157480315):0.0078740157480315,29:0.0078740157480315):0.0078740157480315,30:0.0078740157480315):0.0078740157480315,31:0.0078740157480315):0.0078740157480315,32:0.0078740157480315):0.0078740157480315,33:0.0078740157480315):0.0078740157480315,34:0.0078740157480315):0.0078740157480315,35:0.0078740157480315):0.0078740157480315,36:0.0078740157480315):0.0078740157480315,37:0.0078740157480315):0.0078740157480315,38:0.0078740157480315):0.0078740157480315,39:0.0078740157480315):0.0078740157480315,40:0.0078740157480315):0.0078740157480315,41:0.0078740157480315):0.0078740157480315,42:0.0078740157480315):0.0078740157480315,43:0.0078740157480315):0.0078740157480315,44:0.0078740157480315):0.0078740157480315,45:0.0078740157480315):0.0078740157480315,46:0.0078740157480315):0.0078740157480315,47:0.0078740157480315):0.0078740157480315,48:0.0078740157480315):0.0078740157480315,49:0.0078740157480315):0.0078740157480315,50:0.0078740157480315):0.0078740157480315,51:0.0078740157480315):0.0078740157480315,52:0.0078740157480315):0.0078740157480315,53:0.0078740157480315):0.0078740157480315,54:0.0078740157480315):0.0078740157480315,55:0.0078740157480315):0.0078740157480315,56:0.0078740157480315):0.0078740157480315,57:0.0078740157480315):0.0078740157480315,58:0.0078740157480315):0.0078740157480315,59:0.0078740157480315):0.0078740157480315,60:0.0078740157480315):0.0078740157480315,61:0.0078740157480315):0.0078740157480315,62:0.0078740157480315):0.0078740157480315,63:0.0078740157480315):0.0078740157480315,64:0.0078740157480315):0.0078740157480315,65:0.0078740157480315):0.0078740157480315,66:0.0078740157480315):0.0078740157480315,67:0.0078740157480315):0.0078740157480315,68:0.0078740157480315):0.0078740157480315,69:0.0078740157480315):0.0078740157480315,70:0.0078740157480315):0.0078740157480315,71:0.0078740157480315):0.0078740157480315,72:0.0078740157480315):0.0078740157480315,73:0.0078740157480315):0.0078740157480315,74:0.0078740157480315):0.0078740157480315,75:0.0078740157480315):0.0078740157480315,76:0.0078740157480315):0.0078740157480315,77:0.0078740157480315):0.0078740157480315,78:0.0078740157480315):0.0078740157480315,79:0.0078740157480315):0.0078740157480315,80:0.0078740157480315):0.0078740157480315,81:0.0078740157480315):0.0078740157480315,82:0.0078740157480315):0.0078740157480315,83:0.0078740157480315):0.0078740157480315,84:0.0078740157480315):0.0078740157480315,85:0.0078740157480315):0.0078740157480315,86:0.0078740157480315):0.0078740157480315,87:0.0078740157480315):0.0078740157480315,88:0.0078740157480315):0.0078740157480315,89:0.0078740157480315):0.0078740157480315,90:0.0078740157480315):0.0078740157480315,91:0.0078740157480315):0.0078740157480315,92:0.0078740157480315):0.0078740157480315,93:0.0078740157480315):0.0078740157480315,94:0.0078740157480315):0.0078740157480315,95:0.0078740157480315):0.0078740157480315,96:0.0078740157480315):0.0078740157480315,97:0.0078740157480315):0.0078740157480315,98:0.0078740157480315):0.0078740157480315,99:0.0078740157480315):0.0078740157480315,100:0.0078740157480315):0.0078740157480315,101:0.0078740157480315):0.0078740157480315,102:0.0078740157480315):0.0078740157480315,103:0.0078740157480315):0.0078740157480315,104:0.0078740157480315):0.0078740157480315,105:0.0078740157480315):0.0078740157480315,106:0.0078740157480315):0.0078740157480315,107:0.0078740157480315):0.0078740157480315,108:0.0078740157480315):0.0078740157480315,109:0.0078740157480315):0.0078740157480315,110:0.0078740157480315):0.0078740157480315,111:0.0078740157480315):0.0078740157480315,112:0.0078740157480315):0.0078740157480315,113:0.0078740157480315):0.0078740157480315,114:0.0078740157480315):0.0078740157480315,115:0.0078740157480315):0.0078740157480315,116:0.0078740157480315):0.0078740157480315,117:0.0078740157480315):0.0078740157480315,118:0.0078740157480315):0.0078740157480315,119:0.0078740157480315):0.0078740157480315,120:0.0078740157480315):0.0078740157480315,121:0.0078740157480315):0.0078740157480315,122:0.0078740157480315):0.0078740157480315,123:0.0078740157480315):0.0078740157480315,124:0.0078740157480315):0.0078740157480315,125:0.0078740157480315):0.0078740157480315,126:0.0078740157480315):0.0078740157480315,127:0.0078740157480315):0.0078740157480315,128:0.0078740157480315);"
mp27str <- "(((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((1:0.0078740157480315,2:0.0078740157480315):0.0078740157480315,3:0.0116921783447874):0.0078740157480315,4:0.0157144273508188):0.0078740157480315,5:0.0189070659973398):0.0078740157480315,6:0.0215564536872705):0.0078740157480315,7:0.0238216668561823):0.0078740157480315,8:0.0258005206221736):0.0078740157480315,9:0.0275575564592935):0.0078740157480315,10:0.0291376406755308):0.0078740157480315,11:0.0305732434849321):0.0078740157480315,12:0.0318886299783429):0.0078740157480315,13:0.0331024148588724):0.0078740157480315,14:0.0342291931552556):0.0078740157480315,15:0.0352806214436376):0.0078740157480315,16:0.0362661580023307):0.0078740157480315,17:0.03719358347994):0.0078740157480315,18:0.0380693758823827):0.0078740157480315,19:0.0388989862418931):0.0078740157480315,20:0.0396870449706884):0.0078740157480315,21:0.0404375188274153):0.0078740157480315,22:0.0411538320416336):0.0078740157480315,23:0.0418389609943079):0.0078740157480315,24:0.0424955090961902):0.0078740157480315,25:0.0431257666369267):0.0078740157480315,26:0.0437317590868173):0.0078740157480315,27:0.0443152864266105):0.0078740157480315,28:0.0448779554343439):0.0078740157480315,29:0.0454212063909342):0.0078740157480315,30:0.0459463353239945):0.0078740157480315,31:0.0464545126557864):0.0078740157480315,32:0.04694679893123):0.0078740157480315,33:0.0474241581581312):0.0078740157480315,34:0.0478874691819337):0.0078740157480315,35:0.0483375354326419):0.0078740157480315,36:0.0487750933157507):0.0078740157480315,37:0.0492008194674851):0.0078740157480315,38:0.0496153370539845):0.0078740157480315,39:0.0500192212617605):0.0078740157480315,40:0.050413004100931):0.0078740157480315,41:0.0507971786219501):0.0078740157480315,42:0.0511722026297396):0.0078740157480315,43:0.0515385019654475):0.0078740157480315,44:0.0518964734148707):0.0078740157480315,45:0.0522464872933812):0.0078740157480315,46:0.0525888897495996):0.0078740157480315,47:0.0529240048237635):0.0078740157480315,48:0.0532521362914845):0.0078740157480315,49:0.0535735693192029):0.0078740157480315,50:0.0538885719539565):0.0078740157480315,51:0.0541973964669744):0.0078740157480315,52:0.0545002805679806):0.0078740157480315,53:0.0547974485048535):0.0078740157480315,54:0.0550891120613917):0.0078740157480315,55:0.055375471464312):0.0078740157480315,56:0.0556567162092103):0.0078740157480315,57:0.0559330258140235):0.0078740157480315,58:0.0562045705074973):0.0078740157480315,59:0.0564715118592742):0.0078740157480315,60:0.0567340033574443):0.0078740157480315,61:0.0569921909387276):0.0078740157480315,62:0.057246213475879):0.0078740157480315,63:0.0574962032263875):0.0078740157480315,64:0.0577422862461048):0.0078740157480315,65:0.0579845827710398):0.0078740157480315,66:0.0582232075702157):0.0078740157480315,67:0.0584582702721818):0.0078740157480315,68:0.0586898756675069):0.0078740157480315,69:0.0589181239893435):0.0078740157480315,70:0.0591431111739423):0.0078740157480315,71:0.059364929102814):0.0078740157480315,72:0.0595836658280665):0.0078740157480315,73:0.0597994057823036):0.0078740157480315,74:0.0600122299743345):0.0078740157480315,75:0.0602222161718319):0.0078740157480315,76:0.0604294390719686):0.0078740157480315,77:0.0606339704609707):0.0078740157480315,78:0.0608358793634379):0.0078740157480315,79:0.0610352321822117):0.0078740157480315,80:0.0612320928294986):0.0078740157480315,81:0.0614265228498952):0.0078740157480315,82:0.0616185815359122):0.0078740157480315,83:0.0618083260365347):0.0078740157480315,84:0.06199581145932):0.0078740157480315,85:0.0621810909664868):0.0078740157480315,86:0.062364215865418):0.0078740157480315,87:0.062545235693959):0.0078740157480315,88:0.0627241983008688):0.0078740157480315,89:0.0629011499217512):0.0078740157480315,90:0.0630761352507647):0.0078740157480315,91:0.0632491975083929):0.0078740157480315,92:0.0634203785055298):0.0078740157480315,93:0.0635897187041202):0.0078740157480315,94:0.0637572572745713):0.0078740157480315,95:0.0639230321501447):0.0078740157480315,96:0.0640870800785136):0.0078740157480315,97:0.0642494366706611):0.0078740157480315,98:0.0644101364472848):0.0078740157480315,99:0.0645692128828558):0.0078740157480315,100:0.0647266984474751):0.0078740157480315,101:0.0648826246466577):0.0078740157480315,102:0.0650370220591646):0.0078740157480315,103:0.0651899203730001):0.0078740157480315,104:0.0653413484196763):0.0078740157480315,105:0.0654913342068472):0.0078740157480315,106:0.0656399049494014):0.0078740157480315,107:0.0657870870991033):0.0078740157480315,108:0.0659329063728611):0.0078740157480315,109:0.0660773877796973):0.0078740157480315,110:0.0662205556464941):0.0078740157480315,111:0.0663624336425782):0.0078740157480315,112:0.0665030448032082):0.0078740157480315,113:0.0666424115520223):0.0078740157480315,114:0.0667805557225019):0.0078740157480315,115:0.0669174985785009):0.0078740157480315,116:0.0670532608338907):0.0078740157480315,117:0.0671878626713655):0.0078740157480315,118:0.0673213237604503):0.0078740157480315,119:0.0674536632747526):0.0078740157480315,120:0.0675848999084945):0.0078740157480315,121:0.0677150518923623):0.0078740157480315,122:0.0678441370087062):0.0078740157480315,123:0.0679721726061218):0.0078740157480315,124:0.0680991756134436):0.0078740157480315,125:0.0682251625531793):0.0078740157480315,126:0.06835014955441):0.0078740157480315,127:0.0684741523651827):0.0078740157480315,128:0.0685971863644196);"
mp32str <- mp14str
mp33str <- mp14str
mp34str <- mp14str
newickStrings <- c(sp21str, sp22str, sp23str, sp24str, sp25str, sp26str, sp27str, sp28str, sp29str, sp30str, mp14str, mp15str, mp16str, mp17str, mp18str, mp19str, mp20str, mp21str, mp22str, mp23str, mp24str, mp25str, mp26str, mp27str, mp28str, mp29str, mp30str, mp31str, mp32str, mp33str, mp34str)
names(newickStrings) <- paramSets

sp21bl <- 1/128 # maximally imbalanced, ultrametric
sp22bl <- sp27bl <- sp28bl <- 1/8 # perfectly balanced
sp23bl <- 1/128 # maximally imbalanced, not ultrametric
sp24bl <- 1/128 # maximally imbalanced, not ultrametric
sp25bl <- 1/128 # maximally imbalanced, not ultrametric
sp26bl <- 1/128 # maximally imbalanced, not ultrametric
sp29bl <- sp30bl <- 0 # tiny perfectly balanced ultrametric trees
mp14bl <- sp21bl
mp15bl <- sp22bl
mp16bl <- sp25bl
mp17bl <- sp26bl
mp18bl <- mp14bl
mp19bl <- mp14bl
mp20bl <- mp14bl
mp21bl <- mp18bl
mp22bl <- 1/127
mp23bl <- 1/4
mp24bl <- mp25bl <- mp26bl <- mp27bl <- mp28bl <- mp29bl <- mp30bl <- mp31bl <- 1/127
mp32bl <- mp14bl
mp33bl <- mp14bl
mp34bl <- mp14bl
branchLengths <- c(sp21bl, sp22bl, sp23bl, sp24bl, sp25bl, sp26bl, sp27bl, sp28bl, sp29bl, sp30bl, mp14bl, mp15bl, mp16bl, mp17bl, mp18bl, mp19bl, mp20bl, mp21bl, mp22bl, mp23bl, mp24bl, mp25bl, mp26bl, mp27bl, mp28bl, mp29bl, mp30bl, mp31bl, mp32bl, mp33bl, mp34bl)
names(branchLengths) <- paramSets

dataDir <- "/space/s1/sandra/src/input/treeSim/"
outputDir <- paste0(dataDir, "plots/")
if(!dir.exists(outputDir)) {
  dir.create(outputDir)
}

plotWidth <- 8
plotHeight <- 5.5
sconce2mutColors <- c(SCONCE2="#e98686", SCONCEmut="#22a5e3")

