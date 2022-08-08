library(ggtree)
library(phangorn)
library(ape)
library(reshape2)
library(stringr)
library(ggplot2)
library(cowplot)
library(scales)

source("/space/s1/sandra/src/input/treeSim/scripts/readTreeBranches.R")

# Sat 23 Jul 2022 12:41:50 PM PDT
# same as plotRFdist.R but with mutations
# script to calculate the similarity/distance between estimated trees
# used for comparing nj on t2+t3 distances and some other tree building metric
# plots Robinson-Foulds distances and neighbor joining trees


dataDir <- "/space/s1/sandra/src/input/treeSim/"
outputDir <- paste0(dataDir, "plots/")
if(!dir.exists(outputDir)) {
  dir.create(outputDir)
}
k <- 10
numCellsList <- 20
#numCellsList <- c(20, 40, 60, 80, 100, 120)
filekeys <- c(
              "scAllPMuts_v16", "scAllP_v26")
paramSets <- c("muts/params14", "muts/params15", "muts/params16", "muts/params17", "muts/params32", "muts/params33", "muts/params34")
forceRecalc <- T

treeBranchList <- list()
treeListList <- list()
treeDistList <- list()
njTreesPlotListList <- list()
rfPlotList <- list()

for(numCells in numCellsList) {
  for(paramSet in paramSets) {
    for(key in filekeys) {
      shortName <- paste0(gsub("muts/params", "mp", paramSet), "_", key, "_c", numCells)
      print(paste0("reading ", shortName))
      newickString <- newickStrings[[paramSet]]
      branchLength <- branchLengths[[paramSet]]
      hmmFile <- system(paste0("find ", dataDir, "/", paramSet, " -maxdepth 1 -name \"output_", key, "_", gsub("/", "_", paramSet), "_k", k, "_c", numCells, ".hmm\""), intern=T) # based on scAllP_*sh outBase variable
      if(length(hmmFile) == 0) {
        next
      }
      mergedTreeBranches <- getTreeBranches(newickString, branchLength, hmmFile, forceRecalc)
      if(is.null(mergedTreeBranches)) {
        next
      }
      treeBranchList[[shortName]] <- mergedTreeBranches

      outputFile <- paste0(dataDir, "/plots/njTrees_", gsub("/", "_", paramSet), "_", key, "_k", k, "_c", numCells, "_cnp.png")
      treeList <- createNJtreesCompareMuts(mergedTreeBranches, paramSet, numCells, key)
      if(is.null(treeList)) {
        next
      }
      treeListList[[shortName]] <- treeList

      #njTreesPlotList <- makeNJtreePlots(treeList, outputFile)
      #njTreesPlotListList[[shortName]] <- njTreesPlotList

      dists <- calcRFdistsCompareMuts(treeList)
      # if nearest10, remove t2_t3 from plots (since there might be missing pairs)
      if(grepl("nearest10", key, ignore.case=T)) {
        dists <- subset(dists, variable != "t2_t3")
      }
      dists$paramSet <- paramSet
      treeDistList[[shortName]] <- dists

      p <- makeRFdistPlotMuts(dists, shortName)
      rfPlotList[[shortName]] <- p
    }
  }
}

# euc, cnp2, zzs on sconce, t2+t3
print("plotting just euc/cnp/zzs on sconce, t2+t3")
sconce2MutTreeDistDatList <- list()
sconce2MutTreeDistFiltPlotList <- list()

for(paramSet in paramSets) {
  treeDistDat <- treeDistList[grepl(gsub("muts/params", "mp", paramSet), names(treeDistList))]
  treeDistDat[grepl("Mut", names(treeDistDat))][[1]]$program <- "SCONCEmut"
  treeDistDat[!grepl("Mut", names(treeDistDat))][[1]]$program <- "SCONCE2"
  treeDistDat <- rbind(treeDistDat[[1]], treeDistDat[[2]])
  treeDistDat$program <- factor(treeDistDat$program)
  sconce2MutTreeDistDatList[[paramSet]] <- treeDistDat
  sconce2MutTreeDistFiltPlotList[[paramSet]] <- makeRFdistPlotCompareMuts(treeDistDat, NA)
}

legend <- get_legend(sconce2MutTreeDistFiltPlotList[[1]] + theme(legend.box.margin=margin(0, 0, 0, 0), legend.position="bottom") + guides(colour=guide_legend(nrow=1)))
pGrid <- plot_grid(plotlist=lapply(sconce2MutTreeDistFiltPlotList, FUN=function(x) {x + theme(legend.position="none")}), align='vh', labels="AUTO", ncol=2)
toSave <- plot_grid(pGrid, legend, ncol=1, rel_heights=c(1, 0.05))
outputFile <- paste0(gsub("/", "_", paramSets[1]), "-", gsub("/", "_", paramSets[length(paramSets)]), "_", filekeys[1], "-", filekeys[length(filekeys)], "_k", k, "_c", paste0(numCellsList, collapse="-c"), "_rfDistSconce2SconceMut")
png(paste0(outputDir, "/", outputFile, ".png"), width=plotWidth*1.5, height=plotHeight*2, res=600, units="in"); plot(toSave); dev.off()

