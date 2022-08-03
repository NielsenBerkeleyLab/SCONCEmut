library(ape)
library(plyr)
library(reshape2)
library(stringr)
library(ggplot2)
library(cowplot)

source("/space/s1/sandra/src/input/treeSim/scripts/readTreeBranches.R")

# Mon 13 Jun 2022 12:54:53 PM PDT
# script to make plots of correlations of true branch lengths (from newick trees) and estimated tree branch lengths. Using correlation bc tree branches are scaled differently, but now with mutations!
# runs in ~3 mins

dataDir <- "/space/s1/sandra/src/input/treeSim/"
outputDir <- paste0(dataDir, "plots/")
if(!dir.exists(outputDir)) {
  dir.create(outputDir)
}
k <- 10
numCellsList <- 20
#numCellsList <- c(20,40)
#numCellsList <- c(20, 40, 60, 80, 100, 120)
filekeys <- c(
              #"scAllP_v24_nearest10")
              #"scAllPMuts_v7", "scAllP_v26")
              #"scAllPMuts_v6", "scAllPMuts_v7", "scAllP_v26")
              #"scAllPMuts_v9", "scAllP_v26")
              #"scAllPMuts_v11", "scAllPMuts_v11_singleton", "scAllP_v26")
              #"scAllPMuts_v11", "scAllPMuts_v13", "scAllPMuts_v14", "scAllP_v26")
              #"scAllPMuts_v14", "scAllP_v26")
              #"scAllPMuts_v15", "scAllP_v26")
              "scAllPMuts_v16", "scAllP_v26")
#paramSets <- c("muts/params14", "muts/params15", "muts/params16", "muts/params17")
#paramSets <- c("muts/params16")#, "muts/params17")
#paramSets <- c("muts/params14", "muts/params15", "muts/params16", "muts/params17", "muts/params18", "muts/params19", "muts/params20", "muts/params21")
paramSets <- c("muts/params14", "muts/params15", "muts/params16", "muts/params17", "muts/params32", "muts/params33", "muts/params34")#, "muts/params20")
#paramSets <- c("muts/params14", "muts/params32", "muts/params33", "muts/params34")#, "muts/params20")
#paramSets <- c("muts/params18", "muts/params19", "muts/params20", "muts/params21")
#currMutFilt <- "_anc60_der3"
#mutFilts <- "_mutPairNan"
#mutFilts <- c("", "_mutPairNan", "_mutPairMissing", "_mutPairMissingBranchCap")
mutFilts <- ""
#mutFilts <- paste0("_anc60_der0", c("", paste0("_downSamp", c(500, 1000, 1500, 2000, 5000, 7500, 10000))))
forceRecalc <- T

args <- commandArgs(trailingOnly=TRUE)
if(length(args) > 0) {
  k <- args[1]
  numCellsList <- unlist(strsplit(args[2], ",")) # 20,40,60,80,100,120
  filekeys <- unlist(strsplit(args[3], ",")) # scAllP_v21
  #paramSets <- unlist(strsplit(args[4], ",")) # "segments/params21,segments/params22"
  paramSets <- args[4:length(args)] # "segments/params21" "segments/params22"
}

treeBranchList <- list()
corrPlotList <- list()

for(currMutFilt in mutFilts) {
for(numCells in numCellsList) {
  for(paramSet in paramSets) {
    for(key in filekeys) {
      mutFilt <- ifelse(grepl("Mut", key), currMutFilt, "")
      shortName <- paste0(gsub("muts/params", "mp", paramSet), "_", key, "_c", numCells, mutFilt)
      print(paste0("reading ", shortName))
      newickString <- newickStrings[[paramSet]]
      branchLength <- branchLengths[[paramSet]]
      hmmFile <- system(paste0("find ", dataDir, "/", paramSet, " -name \"output_", key, "_", gsub("/", "_", paramSet), "_k", k, "_c", numCells, mutFilt, ".hmm\""), intern=T) # based on scAllP_*sh outBase variable
      if(grepl("Mut", key) & grepl("params[23]", paramSet)) {
        #hmmFile <- system(paste0("find ", dataDir, "/", paramSet, " -name \"output_", key, "_", "reuseMutEsts_shortcut_", gsub("/", "_", paramSet), "_k", k, "_c", numCells, mutFilt, ".hmm\""), intern=T) # based on scAllP_*sh outBase variable
      }
      if(length(hmmFile) == 0) {
        next
      }
      mergedTreeBranches <- getTreeBranches(newickString, branchLength, hmmFile, forceRecalc)
      if(is.null(mergedTreeBranches)) {
        next
      }
      treeBranchList[[shortName]] <- mergedTreeBranches

      print(paste0("plotting ", shortName))
      p <- makeCorrPlot(mergedTreeBranches, shortName)
      #ggsave(plot=p, file=paste0(outputDir, "/", shortName, "_treeBranchCorr.png"), width=10, height=6) # indv plots
      corrPlotList[[shortName]] <- p
    }
  }
}
}

## separate plot for every cell subset in every paramSet
legend <- get_legend(corrPlotList[[1]] + theme(legend.box.margin=margin(0, 0, 0, 0), legend.position="bottom"))
pGrid <- plot_grid(plotlist=lapply(corrPlotList, FUN=function(x) {x + theme(legend.position="none")}), align='vh', labels="AUTO", nrow=ceiling(length(corrPlotList) / 2), byrow=T)
#pGrid <- plot_grid(plotlist=lapply(corrPlotList, FUN=function(x) {x + theme(legend.position="none")}), align='vh', labels="AUTO", nrow=3, byrow=T)
toSave <- plot_grid(pGrid, legend, ncol=1, rel_heights=c(1, 0.05))
#outputFile <- paste0(gsub("/", "_", paramSets[1]), "-", gsub("/", "_", paramSets[length(paramSets)]), "_", filekeys[1], "-", filekeys[length(filekeys)], "_k", k, "_c", paste0(numCellsList, collapse="-c"), currMutFilt, "_treeBranchCorr")
outputFile <- paste0(gsub("/", "_", paramSets[1]), "-", gsub("/", "_", paramSets[length(paramSets)]), "_", paste0(filekeys, collapse="-"), "_k", k, "_c", paste0(numCellsList, collapse="-c"), currMutFilt, "_treeBranchCorr")
png(paste0(outputDir, "/", outputFile, ".png"), width=4*plotWidth, height=3*plotHeight, res=600, units="in"); plot(toSave); dev.off()

# make sconce2 vs sconceMut treeBranch plots
sconce2MutTreeBranchDatList <- list()
sconce2MutTreeBranchPlotList <- list()
for(paramSet in paramSets) {
  treeBranchDat <- treeBranchList[grepl(gsub("muts/params", "mp", paramSet), names(treeBranchList))]
  treeBranchDat[grepl("Mut", names(treeBranchDat))][[1]]$program <- "SCONCEmut"
  treeBranchDat[!grepl("Mut", names(treeBranchDat))][[1]]$program <- "SCONCE2"
  treeBranchDat <- rbind(treeBranchDat[[1]], treeBranchDat[[2]])
  treeBranchDat$program <- factor(treeBranchDat$program)
  sconce2MutTreeBranchDatList[[paramSet]] <- treeBranchDat
  sconce2MutTreeBranchPlotList[[paramSet]] <- makeCorrPlotCompareMuts(treeBranchDat, NA)
}
legend <- get_legend(corrPlotList[[1]] + theme(legend.box.margin=margin(0, 0, 0, 0), legend.position="bottom"))
pGrid <- plot_grid(plotlist=sconce2MutTreeBranchPlotList[1:4], align='vh', labels="AUTO", ncol=2, byrow=T)
toSave <- plot_grid(pGrid, legend, ncol=1, rel_heights=c(1, 0.05))
outputFile <- paste0(gsub("/", "_", paramSets[1]), "-", gsub("/", "_", paramSets[length(paramSets)]), "_", filekeys[1], "-", filekeys[length(filekeys)], "_k", k, "_c", paste0(numCellsList, collapse="-c"), "_treeBranchSconce2SconceMut_A")
png(paste0(outputDir, "/", outputFile, ".png"), width=plotWidth*1.5, height=plotHeight*1.5, res=600, units="in")
plot(toSave); dev.off()
pGrid <- plot_grid(plotlist=sconce2MutTreeBranchPlotList[5:7], align='vh', labels=c("E", "F", "G"), ncol=2, byrow=T)
toSave <- plot_grid(pGrid, legend, ncol=1, rel_heights=c(1, 0.05))
outputFile <- paste0(gsub("/", "_", paramSets[1]), "-", gsub("/", "_", paramSets[length(paramSets)]), "_", filekeys[1], "-", filekeys[length(filekeys)], "_k", k, "_c", paste0(numCellsList, collapse="-c"), "_treeBranchSconce2SconceMut_B")
png(paste0(outputDir, "/", outputFile, ".png"), width=plotWidth*1.5, height=plotHeight*1.5, res=600, units="in")
plot(toSave); dev.off()



#stop()
#
#
## plot combining cell subsets across each paramSet
#if(length(numCellsList) > 1) {
#  combinedTreeBranchList <- list()
#  combinedCorrPlotList <- list()
#  for(paramSet in paramSets) {
#    for(key in filekeys) {
#      mutFilt <- ifelse(grepl("Mut", key), currMutFilt, "")
#      datasetName <- paste0(paramSet, ", ", key, mutFilt)
#      combinedTreeBranchList[[datasetName]] <- do.call(rbind, treeBranchList[grepl(gsub("muts/params", "mp", paramSet), names(treeBranchList)) & grepl(key, names(treeBranchList))])
#      if(datasetName %in% names(combinedTreeBranchList)) {
#        if(nrow(combinedTreeBranchList[[datasetName]]) != (choose(20,2) * length(numCellsList) * 4)) { # NOTE: skips any incomplete runs. 20 choose 2 pairs * # cell subsets * 4 branch lengths (t1/t2/t3/t2+t3)
#          next
#        }
#        combinedCorrPlotList[[datasetName]] <- makeCorrPlot(combinedTreeBranchList[[datasetName]], datasetName)
#      }
#    }
#  }
#  legend <- get_legend(combinedCorrPlotList[[1]] + theme(legend.box.margin=margin(0, 0, 0, 0), legend.position="bottom"))
#  pGrid <- plot_grid(plotlist=lapply(combinedCorrPlotList, FUN=function(x) {x + theme(legend.position="none")}), align='vh', labels="AUTO", nrow=ceiling(length(combinedCorrPlotList) / 2), byrow=F)
#  toSave <- pGrid
#  #outputFile <- paste0(gsub("/", "_", paramSets[1]), "-", gsub("/", "_", paramSets[length(paramSets)]), "_", filekeys[1], "-", filekeys[length(filekeys)], "_k", k, "_c", paste0(numCellsList, collapse="-c"), mutFilt, "_combinedtreeBranchCorr")
#  outputFile <- paste0(gsub("/", "_", paramSets[1]), "-", gsub("/", "_", paramSets[length(paramSets)]), "_", paste0(filekeys, collapse="-"), "_k", k, "_c", paste0(numCellsList, collapse="-c"), mutFilt, "_combinedtreeBranchCorr")
#  png(paste0(outputDir, "/", outputFile, ".png"), width=2*plotWidth, height=3*plotHeight, res=600, units="in"); plot(toSave); dev.off()
#  #pdf(paste0(outputDir, "/", outputFile, ".pdf"), width=2*plotWidth, height=plotHeight); plot(toSave); dev.off()
#}
#
#
