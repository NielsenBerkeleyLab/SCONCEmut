library(ggplot2)
library(reshape2)
library(cowplot)
library(stringr)

# Fri 18 Feb 2022 06:08:27 PM PST
# file of shared functions for reading bed files from paired analyses

####################################
# SSE
####################################
calcIndvSSE <- function(sconceFileList, meanFileList, medianFileList, modeFileList, groundTruthList, cellIDs) {
  # cell | sconceSumSq | meanSumSq | medianSumSq | modeSumSq
  dat <- do.call(rbind, lapply(cellIDs, FUN=function(cellID) {
    sconceDat <- read.table(sconceFileList[grepl(cellID, sconceFileList, fixed=T)], stringsAsFactors=F, sep="\t", header=F)
    meanDat <- read.table(meanFileList[grepl(cellID, meanFileList, fixed=T)], stringsAsFactors=F, sep="\t", header=F)
    medianDat <- read.table(medianFileList[grepl(cellID, medianFileList, fixed=T)], stringsAsFactors=F, sep="\t", header=F)
    modeDat <- read.table(modeFileList[grepl(cellID, modeFileList, fixed=T)], stringsAsFactors=F, sep="\t", header=F)
    #groundTruth <- read.table(groundTruthFileList[grepl(cellID, groundTruthFileList, fixed=T)], stringsAsFactors=F, sep="\t", header=F)
    groundTruth <- groundTruthList[[cellID]]

    colnames(sconceDat) <- colnames(meanDat) <- colnames(medianDat) <- colnames(modeDat) <- colnames(groundTruth) <- header
    sconceDat$idx <- meanDat$idx <- medianDat$idx <- modeDat$idx <- groundTruth$idx <- 1:nrow(sconceDat)
  
    sconceMerged <- merge(sconceDat, groundTruth, by=c("chr", "start", "end", "idx"))
    meanMerged <- merge(meanDat, groundTruth, by=c("chr", "start", "end", "idx"))
    medianMerged <- merge(medianDat, groundTruth, by=c("chr", "start", "end", "idx"))
    modeMerged <- merge(modeDat, groundTruth, by=c("chr", "start", "end", "idx"))
  
    sconceSumSq <- sum((sconceMerged$copyNumber.x - sconceMerged$copyNumber.y)^2)
    meanSumSq <- sum((meanMerged$copyNumber.x - meanMerged$copyNumber.y)^2)
    medianSumSq <- sum((medianMerged$copyNumber.x - medianMerged$copyNumber.y)^2)
    modeSumSq <- sum((modeMerged$copyNumber.x - modeMerged$copyNumber.y)^2)
    #sconceSumSq <- sum((round(sconceMerged$copyNumber.x) - round(sconceMerged$copyNumber.y))^2)
    #meanSumSq <- sum((round(meanMerged$copyNumber.x) - round(meanMerged$copyNumber.y))^2)
    #medianSumSq <- sum((round(medianMerged$copyNumber.x) - round(medianMerged$copyNumber.y))^2)
    #modeSumSq <- sum((round(modeMerged$copyNumber.x) - round(modeMerged$copyNumber.y))^2)
  
    data.frame(cell=cellID, SCONCE=sconceSumSq, mean=meanSumSq, median=medianSumSq, mode=modeSumSq)
  }))
  dat
}

calcPairsSSE <- function(pairsFileList, groundTruthList, cellIDs) {
  # cells in pair | cell | variable (pair) | value (sumSq)
  dat <- do.call(rbind, lapply(pairsFileList, FUN=function(f) {
    splitFileName <- gsub(".bed", "", gsub(".hg19_lite", "", unlist(str_split(f, "__"))))
    pair <- gsub("pair_", "", splitFileName[2])
    cellID <- str_extract(splitFileName[3], "cancer_cell_[0-9]*.")
    pairDat <- read.table(f, stringsAsFactors=F, sep="\t", header=F)
    colnames(pairDat) <- header
    pairDat$idx <- 1:nrow(pairDat)

    groundTruth <- groundTruthList[[cellID]]
    pairMerged <- merge(pairDat, groundTruth, by=c("chr", "start", "end", "idx"))
    sumSq <- sum((pairMerged$copyNumber.x - pairMerged$copyNumber.y)^2)
    #sumSq <- sum((round(pairMerged$copyNumber.x) - round(pairMerged$copyNumber.y))^2)
    data.frame(pair=pair, cell=cellID, variable="one pair", value=sumSq)
  }))
  dat
}

calcAneuSSE <- function(aneuFileList, groundTruthList, cellIDs, path="/aneu") {
  # cell | variable (aneu) | value (aneuSumSq)
  dat <- do.call(rbind, lapply(cellIDs, FUN=function(cellID) {
    aneuDat <- read.table(aneuFileList[grepl(cellID, aneuFileList, fixed=T)], stringsAsFactors=F, sep="\t", header=F)
    colnames(aneuDat) <- header
    aneuDat$idx <- 1:nrow(aneuDat)
    groundTruth <- groundTruthList[[cellID]]
    aneuMerged <- merge(aneuDat, groundTruth, by=c("chr", "start", "end", "idx"))
    if(nrow(aneuMerged) != nrow(groundTruth)) {
      aneuMerged <- merge(aneuDat, groundTruth, by=c("chr", "start", "idx"))
      if(nrow(aneuMerged) != nrow(groundTruth)) {
        warning(paste0("aneuMerged has the wrong number of rows for ", paramSet, ", ", cellID))
      }
    }
    aneuSumSq <- sum((aneuMerged$copyNumber.x - aneuMerged$copyNumber.y)^2)
    data.frame(cell=cellID, variable="AneuFinder", value=aneuSumSq)
  }))
  dat
}

calcAllSSE <- function(paramSet, key, outputFile, numCells, mutFilt=NULL, forceRecalc=F, inclAneu=F) {
  indvDat <- pairDat <- aneuDat <- NULL
  if(!forceRecalc && file.exists(paste0(outputFile, "_indv.txt")) && file.exists(paste0(outputFile, "_pairs.txt"))) {
    indvDat <- read.table(paste0(outputFile, "_indv.txt"), sep="\t", header=T)
    pairDat <- read.table(paste0(outputFile, "_pairs.txt"), sep="\t", header=T)
    if(inclAneu && file.exists(paste0(outputFile, "_indv.txt"))) {
      aneuDat <- read.table(paste0(outputFile, "_aneu.txt"), sep="\t", header=T)
    }
  } else {
    outputFilePrefix <- paste0("output_", key, "_", gsub("/", "_", paramSet), "_k", k, "_c", numCells, mutFilt) # based on scAllP_*sh outBase variable
    if(grepl("Mut", key) & grepl("params[23]", paramSet)) {
      #outputFilePrefix <- paste0("output_", key, "_", "reuseMutEsts_shortcut_", gsub("/", "_", paramSet), "_k", k, "_c", numCells, mutFilt) # based on scAllP_*sh outBase variable
    }
    unfiltSconceFileList <- system(paste0("find ", dataDir, paramSet, " -maxdepth 1 -name \"", outputFilePrefix, "*__sconce__*.bed\" | sort -V"), intern=T)
    unfiltPairsFileList <- system(paste0("find ", dataDir, paramSet,  " -maxdepth 1 -name \"", outputFilePrefix, "*__pair_*.bed\" | sort -V"), intern=T)
    unfiltMeanFileList <- system(paste0("find ", dataDir, paramSet,   " -maxdepth 1 -name \"", outputFilePrefix, "*__mean.bed\" | sort -V"), intern=T)
    unfiltMedianFileList <- system(paste0("find ", dataDir, paramSet, " -maxdepth 1 -name \"", outputFilePrefix, "*__median.bed\" | sort -V"), intern=T)
    unfiltModeFileList <- system(paste0("find ", dataDir, paramSet,   " -maxdepth 1 -name \"", outputFilePrefix, "*__mode.bed\" | sort -V"), intern=T)
    groundTruthFileList <- system(paste0("find ", dataDir, paramSet,  " -maxdepth 1 -name \"true_cancer_cell_*.bed\" -or -name \"true_copyNumber_cancer_cell_*.bed\" | sort -V"), intern=T)

    # if missing any type of output files to read, return null
    if(length(unfiltSconceFileList) == 0 || length(unfiltPairsFileList) == 0 || length(unfiltMeanFileList) == 0 || length(unfiltMedianFileList) == 0 || length(unfiltModeFileList) == 0) {
      return(NULL)
    }

    tumorDepths <- read.table(paste0(dataDir, paramSet, "/tumor_depths_", numCells), stringsAsFactors=F)$V1
    cellIDs <- str_extract(tumorDepths, "cancer_cell_[0-9]*.")

    sconceFileList <- sapply(cellIDs, FUN=function(cellID) {unfiltSconceFileList[grepl(cellID, unfiltSconceFileList, fixed=T)]})
    pairsFileList <- unfiltPairsFileList[sapply(lapply(str_extract_all(unfiltPairsFileList, "cancer_cell_[0-9]*."), unique), FUN=function(cellNames) {cellNames[1] %in% cellIDs && cellNames[2] %in% cellIDs})]
    meanFileList <- sapply(cellIDs, FUN=function(cellID) {unfiltMeanFileList[grepl(cellID, unfiltMeanFileList, fixed=T)]})
    medianFileList <- sapply(cellIDs, FUN=function(cellID) {unfiltMedianFileList[grepl(cellID, unfiltMedianFileList, fixed=T)]})
    modeFileList <- sapply(cellIDs, FUN=function(cellID) {unfiltModeFileList[grepl(cellID, unfiltModeFileList, fixed=T)]})

    filtGroundTruthFileList <- sapply(cellIDs, FUN=function(cellID) {groundTruthFileList[grepl(cellID, groundTruthFileList, fixed=T)]})
    groundTruthList <- lapply(filtGroundTruthFileList, FUN=function(f) {
      currTruth <- read.table(f, stringsAsFactors=F, sep="\t", header=F)
      colnames(currTruth) <- header
      currTruth$idx <- 1:nrow(currTruth)
      currTruth
    })

    indvDat <- calcIndvSSE(sconceFileList, meanFileList, medianFileList, modeFileList, groundTruthList, cellIDs)
    write.table(indvDat, paste0(outputFile, "_indv.txt"), sep="\t", col.names=T, row.names=F, quote=F)

    pairDat <- calcPairsSSE(pairsFileList, groundTruthList, cellIDs)
    write.table(pairDat, paste0(outputFile, "_pairs.txt"), sep="\t", col.names=T, row.names=F, quote=F)

    if(inclAneu) {
      # find aneu files. assumes $dataDir/$paramSet/aneu/$cell.aneufinderBins.bed
      path <- "/aneu"
      aneuFileList <- sapply(cellIDs, FUN=function(cellID) {
        system(paste0("find ", dataDir, paramSet, path, " -maxdepth 1 -name \"simu_", cellID, "hg19_lite.aneufinderBins.bed\" | sort -V"), intern=T)
      })
      aneuDat <- calcAneuSSE(aneuFileList, groundTruthList, cellIDs, path=path)
      write.table(aneuDat, paste0(outputFile, "_aneu.txt"), sep="\t", col.names=T, row.names=F, quote=F)
    }
  }
  indvDat_m <- melt(indvDat)

  toPlot <- NULL
  if(inclAneu) {
    toPlot <- rbind(indvDat_m, pairDat[,c("cell", "variable", "value")], aneuDat)
    toPlot$variable <- factor(toPlot$variable, levels=c("SCONCE", "one pair", "mean", "median", "mode", "AneuFinder"))
  }
  else {
    toPlot <- rbind(indvDat_m, pairDat[,c("cell", "variable", "value")])
    toPlot$variable <- factor(toPlot$variable, levels=c("SCONCE", "one pair", "mean", "median", "mode"))
  }
  toPlot
}

makeSSEPlot <- function(toPlot, plotTitle) {  
  medians <- do.call(rbind, lapply(levels(toPlot$variable), FUN=function(x) {data.frame(variable=x, median=median(toPlot[toPlot$variable == x, "value"]))}))
  p <- ggplot(toPlot, aes(x=variable, y=value, colour=variable, alpha=(variable == "one pair"))) + geom_boxplot(outlier.alpha=0) + geom_point(position="jitter") + theme_bw() + theme(axis.text.x=element_blank(), legend.title=element_blank()) + stat_summary(fun="median", fun.min="median", fun.max= "median", geom="crossbar", colour="gray35", show.legend=F) + geom_text(data=medians, colour="black", alpha=1, aes(label=sprintf("%.2f", round(median, digits=2)), x=variable),y=Inf, hjust=1, angle=90) + scale_alpha_manual(values=c(1, 0.2), guide="none")
  if(!is.na(plotTitle)) {
    p <- p + labs(x="method", y="SSE", colour="method", title=plotTitle) 
  } else {
    p <- p + labs(x="method", y="SSE", colour="method")
  }
  p
}

makeSSEPlotCompareMuts <- function(toPlot, plotTitle) {  
  medians <- do.call(rbind, lapply(levels(toPlot$variable), FUN=function(x) {
    do.call(rbind, lapply(levels(toPlot$program), FUN=function(prog) {
      data.frame(variable=x, program=prog, median=median(toPlot[toPlot$variable == x & toPlot$program == prog, "value"]))
    }))
  }))
  p <- ggplot(toPlot, aes(x=variable, y=value, colour=program, alpha=(variable == "one pair"))) + geom_boxplot(outlier.alpha=0) + geom_point(position=position_jitterdodge()) + theme_bw() + stat_summary(fun="median", fun.min="median", fun.max= "median", geom="crossbar", colour="gray35", show.legend=F) + geom_text(data=medians, colour="black", alpha=1, aes(label=sprintf("%.2f", round(median, digits=2)), x=variable),y=Inf, hjust=1, vjust=ifelse(medians$program == "SCONCE2", -1, 1), angle=90) + scale_alpha_manual(values=c(1, 0.2), guide="none") + scale_colour_manual(values=sconce2mutColors)
  if(!is.na(plotTitle)) {
    p <- p + labs(x=element_blank(), y="SSE", colour="method", title=plotTitle) 
  } else {
    p <- p + labs(x=element_blank(), y="SSE", colour="method")
  }
  p
}
####################################
# SSE diminishing returns as you add more cells
####################################
getCellOrdering <- function(currSelectionMethod, sconceDistFromCellNum) {
  if(currSelectionMethod == "random") {
    cellSeq <- sample(sconceDistFromCellNum)
  } else if(currSelectionMethod == "nearest") {
    cellSeq <- sort(sconceDistFromCellNum, decreasing=F)
  } else if(currSelectionMethod == "furthest") {
    cellSeq <- sort(sconceDistFromCellNum, decreasing=T)
  }
  cellSeq
}

#calcDiminishingReturns <- function(sconceList, pairsFileList, groundTruthList, btnSconceDist, cellIDs) {
#  # cellID | numCellsSummarized | selectionMethod | sumSq | improvement from sconce | cellAdded
#  dat <- do.call(rbind, lapply(cellIDs, FUN=function(cellID) {
#    sconceDistFromCellNum <- as.matrix(btnSconceDist)[cellID,]
#    sconceDistFromCellNum <- sconceDistFromCellNum[names(sconceDistFromCellNum) != cellID] # remove 0 dist from self
#
#    filtPairsFileList <- pairsFileList[sapply(str_extract_all(pairsFileList, "cancer_cell_[0-9]*."), FUN=function(cellNames) {cellNames[3] == cellID})]
#    pairsList <- lapply(filtPairsFileList, FUN=function(f) {
#      currPair <- read.table(f, stringsAsFactors=F, sep="\t", header=F)
#      colnames(currPair) <- header
#      currPair$idx <- 1:nrow(currPair)
#      currPair
#    })
#    # set the name to be the other cell in the pair (ie that isn't cellID)
#    names(pairsList) <- sapply(str_extract_all(filtPairsFileList, "cancer_cell_[0-9]*."), FUN=function(cellNames) {cellNames[cellNames != cellID]})
#
#    # special case for sconce (ie 1 cell summarized)
#    sconceDat <- sconceList[[cellID]]
#    groundTruth <- groundTruthList[[cellID]]
#    sconceMerged <- merge(sconceDat, groundTruth, by=c("chr", "start", "end", "idx"))
#    sconceSumSq <- sum((sconceMerged$copyNumber.x - sconceMerged$copyNumber.y)^2)
#    sconceEntry <- data.frame(cell=cellID, numCellsSummarized=1, selectionMethod=NA, SSE=sconceSumSq, impr=0, cellAdded=NA)
#
#    # cellID | numCellsSummarized | selectionMethod | sumSq | improvement from sconce | cellAdded
#    sumSqOverCells <- do.call(rbind, lapply(selectionMethods, FUN=function(currSelectionMethod) {
#      currSconceEntry <- sconceEntry
#      currSconceEntry$selectionMethod <- currSelectionMethod
#      cellSeq <- getCellOrdering(currSelectionMethod, sconceDistFromCellNum)
#      pairsCopyNumberMat <- do.call(cbind, lapply(names(cellSeq), FUN=function(selCell) {
#        pairsList[[selCell]]$copyNumber
#      }))
#      colnames(pairsCopyNumberMat) <- names(cellSeq)
#      currSumSqTab <- rbind(currSconceEntry, do.call(rbind, lapply(1:(length(cellIDs)-1), FUN=function(i) {
#        selectedCells <- names(cellSeq[1:i])
#        currPairsCopyNumberMat <- pairsCopyNumberMat[,selectedCells, drop=F]
#        meanCopyNumber <- rowMeans(currPairsCopyNumberMat)
#        sumSq <- sum((groundTruth$copyNumber - meanCopyNumber)^2)
#        cellAdded <- selectedCells[length(selectedCells)]
#        data.frame(cell=cellID, numCellsSummarized=i+1, selectionMethod=currSelectionMethod, SSE=sumSq, impr=sumSq-sconceSumSq, cellAdded=cellAdded) 
#      })))
#      currSumSqTab
#    }))
#    sumSqOverCells
#  }))
#  dat$selectionMethod <- factor(dat$selectionmethod, levels=selectionMethods)
#  dat
#}

calcDiminishingReturns <- function(paramSet, key, outputFile, numCells, forceRecalc=F) {
  changeInSSEDat <- NULL
  if(!forceRecalc && file.exists(paste0(outputFile, ".txt"))) {
    changeInSSEDat <- read.table(paste0(outputFile, ".txt"), sep="\t", header=T)
  } else {
    outputFilePrefix <- paste0("output_", key, "_", gsub("/", "_", paramSet), "_k", k, "_c", numCells) # based on scAllP_*sh outBase variable
    unfiltSconceFileList <- system(paste0("find ", dataDir, paramSet, " -maxdepth 1 -name \"", outputFilePrefix, "__sconce__*.bed\" | sort -V"), intern=T)
    unfiltPairsFileList <- system(paste0("find ", dataDir, paramSet,  " -maxdepth 1 -name \"", outputFilePrefix, "__pair_*.bed\" | sort -V"), intern=T)
    groundTruthFileList <- system(paste0("find ", dataDir, paramSet,  " -maxdepth 1 -name \"true_cancer_cell_*.bed\" | sort -V"), intern=T)

    # if missing any type of output files to read, return null
    if(length(unfiltSconceFileList) == 0 || length(unfiltPairsFileList) == 0) {
      return(NULL)
    }

    tumorDepths <- read.table(paste0(dataDir, paramSet, "/tumor_depths_", numCells), stringsAsFactors=F)$V1
    cellIDs <- str_extract(tumorDepths, "cancer_cell_[0-9]*.")

    sconceFileList <- sapply(cellIDs, FUN=function(cellID) {unfiltSconceFileList[grepl(cellID, unfiltSconceFileList, fixed=T)]})
    pairsFileList <- unfiltPairsFileList[sapply(lapply(str_extract_all(unfiltPairsFileList, "cancer_cell_[0-9]*."), unique), FUN=function(cellNames) {cellNames[1] %in% cellIDs && cellNames[2] %in% cellIDs})]

    filtGroundTruthFileList <- sapply(cellIDs, FUN=function(cellID) {groundTruthFileList[grepl(cellID, groundTruthFileList, fixed=T)]})
    # list of contents of bed files
    groundTruthList <- lapply(filtGroundTruthFileList, FUN=function(f) {
      currTruth <- read.table(f, stringsAsFactors=F, sep="\t", header=F)
      colnames(currTruth) <- header
      currTruth$idx <- 1:nrow(currTruth)
      currTruth
    })

    sconceList <- lapply(sconceFileList, FUN=function(f) {
      currSconce <- read.table(f, stringsAsFactors=F, sep="\t", header=F)
      colnames(currSconce) <- header
      currSconce$idx <- 1:nrow(currSconce)
      currSconce
    })

    # get matrix of sconce calls for measuring distances
    allSconceCallsMat <- do.call(cbind, lapply(sconceList, FUN=function(f) {
      f$copyNumber
    }))
    colnames(allSconceCallsMat) <- str_extract(sconceFileList, "cancer_cell_[0-9]*.")
    btnSconceDist <- dist(t(allSconceCallsMat), method="euclidean")

    #indvDat <- calcIndvSSEAsAddCells(sconceFileList, pairsFileList, groundTruthList, btnSconceDist, cellIDs)
    # cellID | numCellsSummarized | selectionMethod | sumSq | improvement from sconce | cellAdded
    changeInSSEDat <- do.call(rbind, lapply(cellIDs, FUN=function(cellID) {
      sconceDistFromCellNum <- as.matrix(btnSconceDist)[cellID,]
      sconceDistFromCellNum <- sconceDistFromCellNum[names(sconceDistFromCellNum) != cellID] # remove 0 dist from self

      filtPairsFileList <- pairsFileList[sapply(str_extract_all(pairsFileList, "cancer_cell_[0-9]*."), FUN=function(cellNames) {cellNames[3] == cellID})]
      pairsList <- lapply(filtPairsFileList, FUN=function(f) {
        currPair <- read.table(f, stringsAsFactors=F, sep="\t", header=F)
        colnames(currPair) <- header
        currPair$idx <- 1:nrow(currPair)
        currPair
      })
      # set the name to be the other cell in the pair (ie that isn't cellID)
      names(pairsList) <- sapply(str_extract_all(filtPairsFileList, "cancer_cell_[0-9]*."), FUN=function(cellNames) {cellNames[cellNames != cellID]})

      # special case for sconce (ie 1 cell summarized)
      sconceDat <- sconceList[[cellID]]
      groundTruth <- groundTruthList[[cellID]]
      sconceMerged <- merge(sconceDat, groundTruth, by=c("chr", "start", "end", "idx"))
      sconceSumSq <- sum((sconceMerged$copyNumber.x - sconceMerged$copyNumber.y)^2)
      sconceEntry <- data.frame(cell=cellID, numCellsSummarized=1, selectionMethod=NA, SSE=sconceSumSq, impr=0, cellAdded=NA)

      # cellID | numCellsSummarized | selectionMethod | sumSq | improvement from sconce | cellAdded
      sumSqOverCells <- do.call(rbind, lapply(selectionMethods, FUN=function(currSelectionMethod) {
        currSconceEntry <- sconceEntry
        currSconceEntry$selectionMethod <- currSelectionMethod
        cellSeq <- getCellOrdering(currSelectionMethod, sconceDistFromCellNum)
        pairsCopyNumberMat <- do.call(cbind, lapply(names(cellSeq), FUN=function(selCell) {
          pairsList[[selCell]]$copyNumber
        }))
        colnames(pairsCopyNumberMat) <- names(cellSeq)
        currSumSqTab <- rbind(currSconceEntry, do.call(rbind, lapply(1:(length(cellIDs)-1), FUN=function(i) {
          selectedCells <- names(cellSeq[1:i])
          currPairsCopyNumberMat <- pairsCopyNumberMat[,selectedCells, drop=F]
          meanCopyNumber <- rowMeans(currPairsCopyNumberMat)
          sumSq <- sum((groundTruth$copyNumber - meanCopyNumber)^2)
          cellAdded <- selectedCells[length(selectedCells)]
          data.frame(cell=cellID, numCellsSummarized=i+1, selectionMethod=currSelectionMethod, SSE=sumSq, impr=sumSq-sconceSumSq, cellAdded=cellAdded) 
        })))
        currSumSqTab
      }))
      sumSqOverCells
    }))
    changeInSSEDat$selectionMethod <- factor(changeInSSEDat$selectionMethod, levels=selectionMethods)
    write.table(changeInSSEDat, paste0(outputFile, ".txt"), sep="\t", col.names=T, row.names=F, quote=F)
  }
  changeInSSEDat
}

makeDiminishingReturnsPlot <- function(toPlot, plotTitle) {  
  # improvement over sconce, with mean and se crossbars
  p <- ggplot(toPlot, aes(x=numCellsSummarized, y=impr, colour=selectionMethod, group=interaction(cell, selectionMethod))) + theme_bw() + stat_summary(fun.data="mean_se", geom="errorbar", aes(group=selectionMethod), width=0.7, alpha=0.7) + stat_summary(fun=mean, geom="line", aes(group=selectionMethod)) + theme(legend.title=element_blank()) 
  if(!is.na(plotTitle)) {
    p <- p + labs(x=expression(kappa*" = # cells summarized"), y="change in SSE", title=plotTitle)
  } else {
    p <- p + labs(x=expression(kappa*" = # cells summarized"), y="change in SSE")
  }
  p
}

summarizeDiminishingReturns <- function(combinedDimReturnsDat) {
  do.call(rbind, lapply(unique(combinedDimReturnsDat$selectionMethod), FUN=function(selMethod) {
    do.call(rbind, lapply(unique(combinedDimReturnsDat$numCellsSummarized), FUN=function(numCells) {
      currDat <- subset(combinedDimReturnsDat, selectionMethod == selMethod & numCellsSummarized == numCells)
      data.frame(selectionMethod=selMethod, numCellsSummarized=numCells, meanImpr=mean(currDat$impr))
    }))
  }))
}

####################################
# Breakpoints
####################################
calcIndvBreakpoints <- function(sconceFileList, meanFileList, medianFileList, modeFileList, groundTruthList, groundTruthBreakIdxList, numTrueBreakpointsList, cellIDs) {
  # cell | numTrueBreakpoints | variable (source) | num numInferredBreakpoints | summedBreakpointDist
  dat <- do.call(rbind, lapply(cellIDs, FUN=function(cellID) {
    sconceDat <- read.table(sconceFileList[grepl(cellID, sconceFileList, fixed=T)], stringsAsFactors=F, sep="\t", header=F)
    meanDat <- read.table(meanFileList[grepl(cellID, meanFileList, fixed=T)], stringsAsFactors=F, sep="\t", header=F)
    medianDat <- read.table(medianFileList[grepl(cellID, medianFileList, fixed=T)], stringsAsFactors=F, sep="\t", header=F)
    modeDat <- read.table(modeFileList[grepl(cellID, modeFileList, fixed=T)], stringsAsFactors=F, sep="\t", header=F)

    colnames(sconceDat) <- colnames(meanDat) <- colnames(medianDat) <- colnames(modeDat) <- header
    sconceDat$idx <- meanDat$idx <- medianDat$idx <- modeDat$idx <- 1:nrow(sconceDat)
  
    groundTruth <- groundTruthList[[cellID]]
    sconceMerged <- merge(sconceDat, groundTruth, by=c("chr", "start", "end", "idx"))
    meanMerged <- merge(meanDat, groundTruth, by=c("chr", "start", "end", "idx"))
    medianMerged <- merge(medianDat, groundTruth, by=c("chr", "start", "end", "idx"))
    modeMerged <- merge(modeDat, groundTruth, by=c("chr", "start", "end", "idx"))
 
    trueBreakIdx <- groundTruthBreakIdxList[[cellID]]
    numTrueBreakpoints <- numTrueBreakpointsList[[cellID]]
  
    colnames(sconceMerged)[colnames(sconceMerged) == "copyNumber.y"] <- "trueCopyNumber"
    colnames(sconceMerged)[colnames(sconceMerged) == "copyNumber.x"] <- "inferredCopyNumber"
    colnames(meanMerged)[colnames(meanMerged) == "copyNumber.y"] <- "trueCopyNumber"
    colnames(meanMerged)[colnames(meanMerged) == "copyNumber.x"] <- "inferredCopyNumber"
    colnames(medianMerged)[colnames(medianMerged) == "copyNumber.y"] <- "trueCopyNumber"
    colnames(medianMerged)[colnames(medianMerged) == "copyNumber.x"] <- "inferredCopyNumber"
    colnames(modeMerged)[colnames(modeMerged) == "copyNumber.y"] <- "trueCopyNumber"
    colnames(modeMerged)[colnames(modeMerged) == "copyNumber.x"] <- "inferredCopyNumber"
  
    # sort by idx
    sconceMerged <- sconceMerged[with(sconceMerged, order(idx)),]
    meanMerged <- meanMerged[with(meanMerged, order(idx)),]
    medianMerged <- medianMerged[with(medianMerged, order(idx)),]
    modeMerged <- modeMerged[with(modeMerged, order(idx)),]
  
    # calc rle for num breakpoints and dist to nearest breakpoint for each file
    sconceInferredRle <- rle(sconceMerged$inferredCopyNumber)
    sconceNumInferredBreakpoints <- length(sconceInferredRle$values)
    sconceInferredBreakIdx <- head(cumsum(sconceInferredRle$lengths),-1)
    sconceSummedBreakpointDist <- sum(sapply(trueBreakIdx, FUN=function(x) {min(abs(sconceInferredBreakIdx - x))}))
  
    meanInferredRle <- rle(meanMerged$inferredCopyNumber)
    meanNumInferredBreakpoints <- length(meanInferredRle$values)
    meanInferredBreakIdx <- head(cumsum(meanInferredRle$lengths),-1)
    meanSummedBreakpointDist <- sum(sapply(trueBreakIdx, FUN=function(x) {min(abs(meanInferredBreakIdx - x))}))
  
    medianInferredRle <- rle(medianMerged$inferredCopyNumber)
    medianNumInferredBreakpoints <- length(medianInferredRle$values)
    medianInferredBreakIdx <- head(cumsum(medianInferredRle$lengths),-1)
    medianSummedBreakpointDist <- sum(sapply(trueBreakIdx, FUN=function(x) {min(abs(medianInferredBreakIdx - x))}))
  
    modeInferredRle <- rle(modeMerged$inferredCopyNumber)
    modeNumInferredBreakpoints <- length(modeInferredRle$values)
    modeInferredBreakIdx <- head(cumsum(modeInferredRle$lengths),-1)
    modeSummedBreakpointDist <- sum(sapply(trueBreakIdx, FUN=function(x) {min(abs(modeInferredBreakIdx - x))}))
  
    sources <- c("SCONCE", "mean", "median", "mode")
    allNumInferredBreakpoints <- c(sconceNumInferredBreakpoints, meanNumInferredBreakpoints, medianNumInferredBreakpoints, modeNumInferredBreakpoints)
    allSummedBreakpointDist <- c(sconceSummedBreakpointDist, meanSummedBreakpointDist, medianSummedBreakpointDist, modeSummedBreakpointDist)
    data.frame(cell=cellID, numTrueBreakpoints=numTrueBreakpoints, variable=sources, numInferredBreakpoints=allNumInferredBreakpoints, summedBreakpointDist=allSummedBreakpointDist, omega=allNumInferredBreakpoints / numTrueBreakpoints)
  }))
  dat
}

calcPairsBreakpoints <- function(pairsFileList, groundTruthList, groundTruthBreakIdxList, numTrueBreakpointsList, cellIDs) {
  # cells in pair | cell | variable (pair) | numInferredBreakpoints | summedBreakpointDist
  dat <- do.call(rbind, lapply(pairsFileList, FUN=function(f) {
    splitFileName <- gsub(".bed", "", gsub(".hg19_lite", "", unlist(str_split(f, "__"))))
    pair <- gsub("pair_", "", splitFileName[2])
    cellID <- str_extract(splitFileName[3], "cancer_cell_[0-9]*.")
    pairDat <- read.table(f, stringsAsFactors=F, sep="\t", header=F)
    colnames(pairDat) <- header
    pairDat$idx <- 1:nrow(pairDat)

    groundTruth <- groundTruthList[[cellID]]
    pairMerged <- merge(pairDat, groundTruth, by=c("chr", "start", "end", "idx"))
    colnames(pairMerged)[colnames(pairMerged) == "copyNumber.y"] <- "trueCopyNumber"
    colnames(pairMerged)[colnames(pairMerged) == "copyNumber.x"] <- "inferredCopyNumber"

    trueBreakIdx <- groundTruthBreakIdxList[[cellID]]
    numTrueBreakpoints <- numTrueBreakpointsList[[cellID]]

    pairMerged <- pairMerged[with(pairMerged, order(idx)),]
    pairInferredRle <- rle(pairMerged$inferredCopyNumber)
    pairNumInferredBreakpoints <- length(pairInferredRle$values)
    pairInferredBreakIdx <- head(cumsum(pairInferredRle$lengths),-1)
    pairSummedBreakpointDist <- sum(sapply(trueBreakIdx, FUN=function(x) {min(abs(pairInferredBreakIdx - x))}))
 
    data.frame(pair=pair, cell=cellID, numTrueBreakpoints=numTrueBreakpoints, variable="one pair", numInferredBreakpoints=pairNumInferredBreakpoints, summedBreakpointDist=pairSummedBreakpointDist, omega=pairNumInferredBreakpoints / numTrueBreakpoints)
  }))
  dat
}

calcAneuBreakpoints <- function(aneuFileList, groundTruthList, groundTruthBreakIdxList, numTrueBreakpointsList, cellIDs, path="/aneu") {
  # cell | variable (AneuFinder) | numInferredBreakpoints | summedBreakpointDist
  dat <- do.call(rbind, lapply(cellIDs, FUN=function(cellID) {
    aneuDat <- read.table(aneuFileList[grepl(cellID, aneuFileList, fixed=T)], stringsAsFactors=F, sep="\t", header=F)
    colnames(aneuDat) <- header
    aneuDat$idx <- 1:nrow(aneuDat)
    groundTruth <- groundTruthList[[cellID]]
    aneuMerged <- merge(aneuDat, groundTruth, by=c("chr", "start", "end", "idx"))
    if(nrow(aneuMerged) != nrow(groundTruth)) {
      aneuMerged <- merge(aneuDat, groundTruth, by=c("chr", "start", "idx"))
      if(nrow(aneuMerged) != nrow(groundTruth)) {
        warning(paste0("aneuMerged has the wrong number of rows for ", paramSet, ", ", cellID))
      }
    }

    trueBreakIdx <- groundTruthBreakIdxList[[cellID]]
    numTrueBreakpoints <- numTrueBreakpointsList[[cellID]]

    colnames(aneuMerged)[colnames(aneuMerged) == "copyNumber.y"] <- "trueCopyNumber"
    colnames(aneuMerged)[colnames(aneuMerged) == "copyNumber.x"] <- "inferredCopyNumber"

    aneuMerged <- aneuMerged[with(aneuMerged, order(idx)),]
    aneuInferredRle <- rle(aneuMerged$inferredCopyNumber)
    aneuNumInferredBreakpoints <- length(aneuInferredRle$values)
    aneuInferredBreakIdx <- head(cumsum(aneuInferredRle$lengths),-1)
    aneuSummedBreakpointDist <- sum(sapply(trueBreakIdx, FUN=function(x) {min(abs(aneuInferredBreakIdx - x))}))
 
    data.frame(cell=cellID, numTrueBreakpoints=numTrueBreakpoints, variable="AneuFinder", numInferredBreakpoints=aneuNumInferredBreakpoints, summedBreakpointDist=aneuSummedBreakpointDist, omega=aneuNumInferredBreakpoints / numTrueBreakpoints)
  }))
  dat
}

calcAllBreakpoints <- function(paramSet, key, outputFile, numCells, mutFilt=NULL, forceRecalc=F, inclAneu=F) {
  indvDat <- pairDat <- aneuDat <- NULL
  if(!forceRecalc && file.exists(paste0(outputFile, "_indv.txt")) && file.exists(paste0(outputFile, "_pairs.txt"))) {
    indvDat <- read.table(paste0(outputFile, "_indv.txt"), sep="\t", header=T)
    pairDat <- read.table(paste0(outputFile, "_pairs.txt"), sep="\t", header=T)
    if(inclAneu && file.exists(paste0(outputFile, "_indv.txt"))) {
      aneuDat <- read.table(paste0(outputFile, "_aneu.txt"), sep="\t", header=T)
    }
  } else {
    outputFilePrefix <- paste0("output_", key, "_", gsub("/", "_", paramSet), "_k", k, "_c", numCells, mutFilt) # based on scAllP_*sh outBase variable
    if(grepl("Mut", key) & grepl("params[23]", paramSet)) {
      #outputFilePrefix <- paste0("output_", key, "_", "reuseMutEsts_shortcut_", gsub("/", "_", paramSet), "_k", k, "_c", numCells, mutFilt) # based on scAllP_*sh outBase variable
    }
    unfiltSconceFileList <- system(paste0("find ", dataDir, paramSet, " -maxdepth 1 -name \"", outputFilePrefix, "*__sconce__*.bed\" | sort -V"), intern=T)
    unfiltPairsFileList <- system(paste0("find ", dataDir, paramSet,  " -maxdepth 1 -name \"", outputFilePrefix, "*__pair_*.bed\" | sort -V"), intern=T)
    unfiltMeanFileList <- system(paste0("find ", dataDir, paramSet,   " -maxdepth 1 -name \"", outputFilePrefix, "*__mean.bed\" | sort -V"), intern=T)
    unfiltMedianFileList <- system(paste0("find ", dataDir, paramSet, " -maxdepth 1 -name \"", outputFilePrefix, "*__median.bed\" | sort -V"), intern=T)
    unfiltModeFileList <- system(paste0("find ", dataDir, paramSet,   " -maxdepth 1 -name \"", outputFilePrefix, "*__mode.bed\" | sort -V"), intern=T)
    groundTruthFileList <- system(paste0("find ", dataDir, paramSet,  " -maxdepth 1 -name \"true_cancer_cell_*.bed\" -or -name \"true_copyNumber_cancer_cell_*.bed\" | sort -V"), intern=T)

    # if missing any type of output files to read, return null
    if(length(unfiltSconceFileList) == 0 || length(unfiltPairsFileList) == 0 || length(unfiltMeanFileList) == 0 || length(unfiltMedianFileList) == 0 || length(unfiltModeFileList) == 0) {
      return(NULL)
    }

    tumorDepths <- read.table(paste0(dataDir, paramSet, "/tumor_depths_", numCells), stringsAsFactors=F)$V1
    cellIDs <- str_extract(tumorDepths, "cancer_cell_[0-9]*.")

    sconceFileList <- sapply(cellIDs, FUN=function(cellID) {unfiltSconceFileList[grepl(cellID, unfiltSconceFileList, fixed=T)]})
    pairsFileList <- unfiltPairsFileList[sapply(lapply(str_extract_all(unfiltPairsFileList, "cancer_cell_[0-9]*."), unique), FUN=function(cellNames) {cellNames[1] %in% cellIDs && cellNames[2] %in% cellIDs})]
    meanFileList <- sapply(cellIDs, FUN=function(cellID) {unfiltMeanFileList[grepl(cellID, unfiltMeanFileList, fixed=T)]})
    medianFileList <- sapply(cellIDs, FUN=function(cellID) {unfiltMedianFileList[grepl(cellID, unfiltMedianFileList, fixed=T)]})
    modeFileList <- sapply(cellIDs, FUN=function(cellID) {unfiltModeFileList[grepl(cellID, unfiltModeFileList, fixed=T)]})

    filtGroundTruthFileList <- sapply(cellIDs, FUN=function(cellID) {groundTruthFileList[grepl(cellID, groundTruthFileList, fixed=T)]})
    groundTruthList <- lapply(filtGroundTruthFileList, FUN=function(f) {
      currTruth <- read.table(f, stringsAsFactors=F, sep="\t", header=F)
      colnames(currTruth) <- header
      currTruth$idx <- 1:nrow(currTruth)
      currTruth
    })
    groundTruthBreakIdxList <- list()
    numTrueBreakpointsList <- list()
    for(cellID in cellIDs) {
      groundTruth <- groundTruthList[[cellID]]
      groundTruthRle <- rle(groundTruth$copyNumber)
      trueBreakIdx <- head(cumsum(groundTruthRle$lengths),-1)
      # if no CNAs, then there aren't any breakpoints
      if(length(trueBreakIdx) == 0) {
        trueBreakIdx <- cumsum(groundTruthRle$lengths)
      }
      numTrueBreakpoints <- length(groundTruthRle$values)
      groundTruthBreakIdxList[[cellID]] <- trueBreakIdx
      numTrueBreakpointsList[[cellID]] <- numTrueBreakpoints
    }

    indvDat <- calcIndvBreakpoints(sconceFileList, meanFileList, medianFileList, modeFileList, groundTruthList, groundTruthBreakIdxList, numTrueBreakpointsList, cellIDs)
    write.table(indvDat, paste0(outputFile, "_indv.txt"), sep="\t", col.names=T, row.names=F, quote=F)

    pairDat <- calcPairsBreakpoints(pairsFileList, groundTruthList, groundTruthBreakIdxList, numTrueBreakpointsList, cellIDs)
    write.table(pairDat, paste0(outputFile, "_pairs.txt"), sep="\t", col.names=T, row.names=F, quote=F)

    if(inclAneu) {
      # find aneu files. assumes $dataDir/$paramSet/aneu/$cell.aneufinderBins.bed
      path <- "/aneu"
      aneuFileList <- sapply(cellIDs, FUN=function(cellID) {
        system(paste0("find ", dataDir, paramSet, path, " -maxdepth 1 -name \"simu_", cellID, "hg19_lite.aneufinderBins.bed\" | sort -V"), intern=T)
      })
      aneuDat <- calcAneuBreakpoints(aneuFileList, groundTruthList, groundTruthBreakIdxList, numTrueBreakpointsList, cellIDs, path=path)
      write.table(aneuDat, paste0(outputFile, "_aneu.txt"), sep="\t", col.names=T, row.names=F, quote=F)
    }
  }

  toPlot <- NULL
  if(inclAneu) {
    toPlot <- rbind(indvDat, pairDat[,-1], aneuDat)
    toPlot$variable <- factor(toPlot$variable, levels=c("SCONCE", "one pair", "mean", "median", "mode", "AneuFinder"))
  }
  else {
    toPlot <- rbind(indvDat, pairDat[,-1])
    toPlot$variable <- factor(toPlot$variable, levels=c("SCONCE", "one pair", "mean", "median", "mode"))
  }
  toPlot
}

makeBreakpointPlot <- function(toPlot, plotTitle) {
  p <- ggplot(toPlot, aes(x=numInferredBreakpoints/numTrueBreakpoints, y=summedBreakpointDist, colour=variable, alpha=(variable == "one pair"))) + geom_point() + theme_bw() + guides(alpha="none") + scale_alpha_manual(values=c(1, 0.1), guide="none") + scale_y_log10() + theme(legend.title=element_blank()) + geom_vline(xintercept=1, linetype="dashed", colour="red", alpha=0.5)
  if(!is.na(plotTitle)) {
    p <- p + labs(x=expression(omega), y="breakpoint dist", title=plotTitle) 
  } else {
    p <- p + labs(x=expression(omega), y="breakpoint dist")
  }
  p
}
makeBreakpointPlotCompareMuts <- function(toPlot, plotTitle) {
  #p <- ggplot(toPlot, aes(x=numInferredBreakpoints/numTrueBreakpoints, y=summedBreakpointDist, colour=variable, shape=program, alpha=(variable == "one pair"))) + geom_point(size=4) + theme_bw() + guides(alpha="none") + scale_alpha_manual(values=c(1, 0.1), guide="none") + scale_y_log10() + theme(legend.title=element_blank()) + geom_vline(xintercept=1, linetype="dashed", colour="red", alpha=0.5) + geom_point(colour="white", size=1.5) + geom_line(aes(group=interaction(cell, variable)))
  p <- ggplot(toPlot, aes(x=numInferredBreakpoints/numTrueBreakpoints, y=summedBreakpointDist, colour=variable, shape=program, alpha=(variable == "one pair"))) + geom_point() + theme_bw() + guides(alpha="none") + scale_alpha_manual(values=c(1, 0.1), guide="none") + scale_y_log10() + theme(legend.title=element_blank()) + geom_vline(xintercept=1, linetype="dashed", colour="red", alpha=0.5) + geom_line(aes(group=interaction(cell, variable)))
  if(!is.na(plotTitle)) {
    p <- p + labs(x=expression(omega), y="breakpoint dist", title=plotTitle) 
  } else {
    p <- p + labs(x=expression(omega), y="breakpoint dist")
  }
  p
}

# write median breakpoint distance and omega values to text and tex files for easy copy/paste into latex tables
writeMedianDistOmegaTexFiles <- function(breakpointDat, medianDistOmegaTexFile) {
  programs <- levels(breakpointDat$variable)
  currParamSets <- unique(breakpointDat$paramSet)

  medianDists <- as.data.frame(t(sapply(currParamSets, FUN=function(set) {
    sapply(programs, FUN=function(program) {
      median(subset(breakpointDat, paramSet == set & variable == program)$summedBreakpointDist)
    })
  })))

  medianOmegas <- as.data.frame(t(sapply(currParamSets, FUN=function(set) {
    sapply(programs, FUN=function(program) {
      median(subset(breakpointDat, paramSet == set & variable == program)$omega)
    })
  })))

  write.table(medianDists, file=paste0(medianDistOmegaTexFile, "_medianDists.txt"), sep="\t", quote=F, row.names=T, col.names=NA)
  write.table(medianOmegas, file=paste0(medianDistOmegaTexFile, "_medianOmegas.txt"), sep="\t", quote=F, row.names=T, col.names=NA)

  #texColnames <- c("\\textbf{SCONCE}", "\\textbf{one pair}", "\\textbf{mean}", "\\textbf{median}", "\\textbf{mode}", "\\textbf{AneuFinder}") # may need to insert line breaks so the tables fit nicely. see plotJointBreakpointComparison.R for syntax
  texColnames <- paste0("\\textbf{", programs, "}")
  names(texColnames) <- levels(programs)

  medianDistsTex <- rbind(texColnames, round(medianDists, digits=4))
  medianDistsTex$paramSet <- paste0(" \\\\ \\hline % ", c("", rownames(medianDists)))
  medianDistsTex <- cbind(c("", LETTERS[1:nrow(medianDists)]), medianDistsTex)
  write.table(paste0(apply(medianDistsTex[,-ncol(medianDistsTex)], 1, FUN=function(x) {paste0(x, collapse=" & ")}),  medianDistsTex[,ncol(medianDistsTex)]), file=paste0(medianDistOmegaTexFile, "_medianDists.tex"), col.names=F, row.names=F, quote=F)

  medianOmegasTex <- rbind(texColnames, round(medianOmegas, digits=4))
  medianOmegasTex$paramSet <- paste0(" \\\\ \\hline % ", c("", rownames(medianOmegas)))
  medianOmegasTex <- cbind(c("", LETTERS[1:nrow(medianOmegas)]), medianOmegasTex)
  write.table(paste0(apply(medianOmegasTex[,-ncol(medianOmegasTex)], 1, FUN=function(x) {paste0(x, collapse=" & ")}),  medianOmegasTex[,ncol(medianOmegasTex)]), file=paste0(medianDistOmegaTexFile, "_medianOmegas.tex"), col.names=F, row.names=F, quote=F)
}


####################################
# shared constants
####################################
header <- c("chr", "start", "end", "copyNumber")
dataDir <- "/space/s1/sandra/src/input/treeSim/"
outputDir <- paste0(dataDir, "plots/")
if(!dir.exists(outputDir)) {
  dir.create(outputDir)
}

plotWidth <- 8
plotHeight <- 5.5
sconce2mutColors <- c(SCONCE2="#e98686", SCONCEmut="#22a5e3")

