
# Fri 17 Jun 2022 12:16:27 PM PDT
# script to plot the number of true/simulated mutations vs number of observed mutations (based on read depth) vs number of inferred mutations
# for both individual cells and for pairs

source("/space/s1/sandra/src/input/treeSim/scripts/readTreeBranches.R")
dataDir <- "/space/s1/sandra/src/input/treeSim/"
outputDir <- paste0(dataDir, "plots/")
if(!dir.exists(outputDir)) {
  dir.create(outputDir)
}
k <- 10
numCellsList <- 20
#numCellsList <- c(20, 40, 60, 80, 100, 120)
filekeys <- c(
              "scAllPMuts_v16")
paramSets <- c("muts/params14", "muts/params15", "muts/params16", "muts/params17", "muts/params32", "muts/params33", "muts/params34")
mutFilts <- ""
forceRecalc <- T

inferredIndvMutCountsList <- list()
inferredPairedMutCountsList <- list()

trueMutCountsList <- list()
obsMutCountsList <- list()
corrPlotList <- list()
numMutsVsBranchesIndvPlotList <- list()
numMutsVsBranchesPairedPlotList <- list()
numMutsInferredVsObsVsTrueIndvPlotList <- list()
numMutsInferredVsObsVsTruePairedPlotList <- list()
numMutsInferredVsObsPairedPlotList <- list()


for(mutFilt in mutFilts) {
  for(paramSet in paramSets) {
    for(key in filekeys) {
      datasetName <- paste0(gsub("muts/params", "mp", paramSet), "_", key, mutFilt)
      print(paste0("reading ", datasetName))
      newickString <- newickStrings[[paramSet]]
      branchLength <- branchLengths[[paramSet]]

      currIndvMutCounts <- list()
      currPairedMutCounts <- list()
      for(numCells in numCellsList) {
        print("getting inferred mut counts")
        hmmFile <- system(paste0("find ", dataDir, "/", paramSet, " -name \"output_", key, "_", gsub("/", "_", paramSet), "_k", k, "_c", numCells, mutFilt, ".hmm\""), intern=T) # based on scAllP_*sh outBase variable
        if(length(hmmFile) == 0) {
          next
        }
        inferredMutCounts <- getInferredPerBranchMutCounts(newickString, branchLength, hmmFile, forceRecalc)
        if(is.null(inferredMutCounts)) {
          next
        }
        currIndvMutCounts[[numCells]] <- inferredMutCounts[["indv"]]
        currPairedMutCounts[[numCells]] <- inferredMutCounts[["paired"]]
      }
      inferredIndvMutCountsList[[datasetName]] <- do.call(rbind, currIndvMutCounts)
      inferredPairedMutCountsList[[datasetName]] <- do.call(rbind, currPairedMutCounts)
      if(length(inferredIndvMutCountsList[[datasetName]]) == 0 && length(inferredPairedMutCountsList[[datasetName]]) == 0) {
        next
      }

      print("getting true mut counts")
      trueMutCounts <- getTruePerBranchMutCounts(newickString, branchLength, paramSet, forceRecalc=forceRecalc, uniqOnly=F, numCellsList=numCellsList)
      if(is.null(trueMutCounts)) {
        next
      }
      trueMutCountsList[[datasetName]] <- trueMutCounts

      print("getting obs mut counts")
      obsMutCounts <- getObservedPerBranchMutCounts(newickString, branchLength, paramSet, "", forceRecalc=forceRecalc, uniqOnly=F, numCellsList=numCellsList)
      if(is.null(obsMutCounts)) {
        next
      }
      obsMutCountsList[[datasetName]] <- obsMutCounts
      # save observed mut counts to files
      saveMutCountsToFile(paramSet, numCells, trueMutCounts, obsMutCounts)

      trueMutCounts <- trueMutCountsList[[datasetName]]
      obsMutCounts <- obsMutCountsList[[datasetName]]
      print(paste0("plotting ", datasetName))
      #indvCorrPlots <- makeIndvTrueObsInferredCorrPlot(trueMutCounts[["indv"]], obsMutCounts[["indv"]], inferredIndvMutCountsList[[datasetName]], datasetName)
      pairedCorrPlots <- makePairedTrueObsInferredCorrPlot(trueMutCounts[["paired"]], obsMutCounts[["paired"]], inferredPairedMutCountsList[[datasetName]], NULL)
      #corrPlotList[[datasetName]] <- list(indvCorrPlots=indvCorrPlots, pairedCorrPlots=pairedCorrPlots)
      corrPlotList[[datasetName]] <- list(pairedCorrPlots=pairedCorrPlots)

      #numMutsVsBranchesIndvPlotList[[datasetName]] <- indvCorrPlots[["pNumMutsVsBranches"]]
      numMutsVsBranchesPairedPlotList[[datasetName]] <- pairedCorrPlots[["pNumMutsVsBranches"]]
      #numMutsInferredVsObsVsTrueIndvPlotList[[datasetName]] <-indvCorrPlots[["mutCorr"]]
      numMutsInferredVsObsVsTruePairedPlotList[[datasetName]] <- pairedCorrPlots[["mutCorr"]]
      numMutsInferredVsObsPairedPlotList[[datasetName]] <- makePairedObsInferredCorrPlot(obsMutCounts[["paired"]], inferredPairedMutCountsList[[datasetName]], NULL)
    }
  }
}

# suppl figure, broken into 2 figures, 4 and 3 paramSets each
legend <- get_legend(numMutsVsBranchesPairedPlotList[[1]] + theme(legend.box.margin=margin(0, 0, 0, 0), legend.position="bottom"))
pGrid <- plot_grid(plotlist=lapply(numMutsVsBranchesPairedPlotList[1:4], FUN=function(x) {x + theme(legend.position="none") + facet_grid(variable ~ branch, scales="free", labeller=label_parsed)}), align='vh', labels="AUTO", ncol=2, byrow=T)
toSave <- plot_grid(pGrid, legend, ncol=1, rel_heights=c(1, 0.05))
outputFile <- paste0(gsub("/", "_", paramSets[1]), "-", gsub("/", "_", paramSets[length(paramSets)]), "_", filekeys[1], "-", filekeys[length(filekeys)], "_k", k, "_c", paste0(numCellsList, collapse="-c"), mutFilt, "_numMutsVsBranchesPaired_A")
png(paste0(outputDir, "/", outputFile, ".png"), width=plotWidth*1.5, height=plotHeight*2, res=600, units="in"); plot(toSave); dev.off()
pGrid <- plot_grid(plotlist=lapply(numMutsVsBranchesPairedPlotList[5:7], FUN=function(x) {x + theme(legend.position="none") + facet_grid(variable ~ branch, scales="free", labeller=label_parsed)}), align='vh', labels=c("E", "F", "G"), ncol=2, byrow=T)
toSave <- plot_grid(pGrid, legend, ncol=1, rel_heights=c(1, 0.05))
outputFile <- paste0(gsub("/", "_", paramSets[1]), "-", gsub("/", "_", paramSets[length(paramSets)]), "_", filekeys[1], "-", filekeys[length(filekeys)], "_k", k, "_c", paste0(numCellsList, collapse="-c"), mutFilt, "_numMutsVsBranchesPaired_B")
png(paste0(outputDir, "/", outputFile, ".png"), width=plotWidth*1.5, height=plotHeight*2*3/4, res=600, units="in"); plot(toSave); dev.off()
#png(paste0(outputDir, "/", outputFile, ".png"), width=25, height=20, res=600, units="in"); plot(toSave); dev.off()

# suppl figure, broken into 2 figures, 4 and 3 paramSets each
pGrid <- plot_grid(plotlist=lapply(numMutsInferredVsObsVsTruePairedPlotList[1:4], FUN=function(x) {x + theme(legend.position="none")}), align='vh', labels="AUTO", ncol=1, byrow=T)
toSave <- plot_grid(pGrid, legend, ncol=1, rel_heights=c(1, 0.05))
outputFile <- paste0(gsub("/", "_", paramSets[1]), "-", gsub("/", "_", paramSets[length(paramSets)]), "_", filekeys[1], "-", filekeys[length(filekeys)], "_k", k, "_c", paste0(numCellsList, collapse="-c"), mutFilt, "_numMutsInferredVsObsVsTruePaired_A")
png(paste0(outputDir, "/", outputFile, ".png"), width=plotWidth*1.75, height=plotHeight*2, res=600, units="in"); plot(toSave); dev.off()
pGrid <- plot_grid(plotlist=lapply(numMutsInferredVsObsVsTruePairedPlotList[5:7], FUN=function(x) {x + theme(legend.position="none")}), align='vh', labels=c("E", "F", "G"), ncol=1, byrow=T)
toSave <- plot_grid(pGrid, legend, ncol=1, rel_heights=c(1, 0.05))
outputFile <- paste0(gsub("/", "_", paramSets[1]), "-", gsub("/", "_", paramSets[length(paramSets)]), "_", filekeys[1], "-", filekeys[length(filekeys)], "_k", k, "_c", paste0(numCellsList, collapse="-c"), mutFilt, "_numMutsInferredVsObsVsTruePaired_B")
png(paste0(outputDir, "/", outputFile, ".png"), width=plotWidth*1.75, height=plotHeight*2*3/4, res=600, units="in"); plot(toSave); dev.off()

# pull out inferred vs obs muts for main figure
legend <- get_legend(numMutsInferredVsObsPairedPlotList[[1]] + theme(legend.box.margin=margin(0, 0, 0, 0), legend.position="bottom"))
pGrid <- plot_grid(plotlist=lapply(numMutsInferredVsObsPairedPlotList, FUN=function(x) {x + theme(legend.position="none")}), align='vh', labels="AUTO", ncol=2, byrow=T)
toSave <- plot_grid(pGrid, legend, ncol=1, rel_heights=c(1, 0.05))
outputFile <- paste0(gsub("/", "_", paramSets[1]), "-", gsub("/", "_", paramSets[length(paramSets)]), "_", filekeys[1], "-", filekeys[length(filekeys)], "_k", k, "_c", paste0(numCellsList, collapse="-c"), mutFilt, "_numMutsInferredVsObsPaired")
png(paste0(outputDir, "/", outputFile, ".png"), width=plotWidth*1.5, height=plotHeight*2, res=600, units="in"); plot(toSave); dev.off()

