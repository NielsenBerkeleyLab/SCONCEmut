
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
              #"scAllP_v24_nearest10")
              #"scAllPMuts_v7", "scAllP_v26")
              #"scAllPMuts_v6", "scAllPMuts_v7")
              #"scAllPMuts_v11", "scAllPMuts_v11_singleton")
              #"scAllPMuts_v11")
              #"scAllPMuts_v14")
              #"scAllPMuts_v15")
              "scAllPMuts_v16")
#paramSets <- c("muts/params14", "muts/params15", "muts/params16", "muts/params17")
#paramSets <- c("muts/params14", "muts/params16", "muts/params17")
#paramSets <- c("muts/params16")
#paramSets <- c("muts/params18", "muts/params19")
#paramSets <- c("muts/params14", "muts/params15", "muts/params16", "muts/params17", "muts/params18", "muts/params19", "muts/params20", "muts/params21")
paramSets <- c("muts/params14", "muts/params15", "muts/params16", "muts/params17", "muts/params32", "muts/params33", "muts/params34") #, "muts/params20")
#paramSets <- c("muts/params14", "muts/params32", "muts/params33", "muts/params34")#, "muts/params20")
#paramSets <- c("muts/params33")
#mutFilts <- c(".anc3_der1", ".anc1_der1", ".anc2_der1", ".anc60_der3", ".full")
#mutFilts <- c("", "_mutPairNan", "_mutPairMissing", "_mutPairMissingBranchCap")
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
    #mutFilt <- ifelse(grepl("Mut", key), "_anc60_der3", "")
    #mutFilt <- "_mutPairNan"
    datasetName <- paste0(gsub("muts/params", "mp", paramSet), "_", key, mutFilt)
    print(paste0("reading ", datasetName))
    newickString <- newickStrings[[paramSet]]
    branchLength <- branchLengths[[paramSet]]

    currIndvMutCounts <- list()
    currPairedMutCounts <- list()
    for(numCells in numCellsList) {
      print("getting inferred mut counts")
      hmmFile <- system(paste0("find ", dataDir, "/", paramSet, " -name \"output_", key, "_", gsub("/", "_", paramSet), "_k", k, "_c", numCells, mutFilt, ".hmm\""), intern=T) # based on scAllP_*sh outBase variable
      if(grepl("Mut", key) & grepl("params[23]", paramSet)) {
        #hmmFile <- system(paste0("find ", dataDir, "/", paramSet, " -name \"output_", key, "_", "reuseMutEsts_shortcut_", gsub("/", "_", paramSet), "_k", k, "_c", numCells, mutFilt, ".hmm\""), intern=T) # based on scAllP_*sh outBase variable
      }
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
    #if(length(currIndvMutCounts) != max(numCellsList) && length(currPairedMutCounts) != max(numCellsList)) { # NOTE: skips any incomplete runs
    #  next
    #}

    print("getting true mut counts")
    trueMutCounts <- getTruePerBranchMutCounts(newickString, branchLength, paramSet, forceRecalc=forceRecalc, uniqOnly=F, numCellsList=numCellsList)
    if(is.null(trueMutCounts)) {
      next
    }
    trueMutCountsList[[datasetName]] <- trueMutCounts

    print("getting obs mut counts")
    #obsMutCounts <- getObservedPerBranchMutCounts(newickString, branchLength, paramSet, mutFilt, forceRecalc=forceRecalc, uniqOnly=F, numCellsList=numCellsList)
    obsMutCounts <- getObservedPerBranchMutCounts(newickString, branchLength, paramSet, "", forceRecalc=forceRecalc, uniqOnly=F, numCellsList=numCellsList)
    if(is.null(obsMutCounts)) {
      next
    }
    obsMutCountsList[[datasetName]] <- obsMutCounts
    # save observed mut counts to files
    saveMutCountsToFile(paramSet, numCells, trueMutCounts, obsMutCounts)
    #stop()

    #next
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

## save plots across datasets
#legend <- get_legend(numMutsVsBranchesIndvPlotList[[1]] + theme(legend.box.margin=margin(0, 0, 0, 0), legend.position="bottom"))
#pGrid <- plot_grid(plotlist=lapply(numMutsVsBranchesIndvPlotList, FUN=function(x) {x + theme(legend.position="none")}), align='vh', labels="AUTO", nrow=ceiling(length(numMutsVsBranchesIndvPlotList) / 2), byrow=F)
#toSave <- plot_grid(pGrid, legend, ncol=1, rel_heights=c(1, 0.05))
#outputFile <- paste0(gsub("/", "_", paramSets[1]), "-", gsub("/", "_", paramSets[length(paramSets)]), "_", filekeys[1], "-", filekeys[length(filekeys)], "_k", k, "_c", paste0(numCellsList, collapse="-c"), mutFilt, "_numMutsVsBranchesIndv")
#png(paste0(outputDir, "/", outputFile, ".png"), width=25, height=20, res=600, units="in"); plot(toSave); dev.off()

#legend <- get_legend(numMutsInferredVsObsVsTrueIndvPlotList[[1]] + theme(legend.box.margin=margin(0, 0, 0, 0), legend.position="bottom"))
#pGrid <- plot_grid(plotlist=lapply(numMutsInferredVsObsVsTrueIndvPlotList, FUN=function(x) {x + theme(legend.position="none")}), align='vh', labels="AUTO", nrow=ceiling(length(numMutsInferredVsObsVsTrueIndvPlotList) / 2), byrow=F)
#toSave <- plot_grid(pGrid, legend, ncol=1, rel_heights=c(1, 0.05))
#outputFile <- paste0(gsub("/", "_", paramSets[1]), "-", gsub("/", "_", paramSets[length(paramSets)]), "_", filekeys[1], "-", filekeys[length(filekeys)], "_k", k, "_c", paste0(numCellsList, collapse="-c"), mutFilt, "_numMutsInferredVsObsVsTrueIndv")
#png(paste0(outputDir, "/", outputFile, ".png"), width=25, height=20, res=600, units="in"); plot(toSave); dev.off()


# suppl figure, broken into 2 figures, 4 and 3 paramSets each
legend <- get_legend(numMutsVsBranchesPairedPlotList[[1]] + theme(legend.box.margin=margin(0, 0, 0, 0), legend.position="bottom"))
#pGrid <- plot_grid(plotlist=lapply(numMutsVsBranchesPairedPlotList, FUN=function(x) {x + theme(legend.position="none") + facet_grid(variable ~ branch, scales="free")}), align='vh', labels="AUTO", nrow=ceiling(length(numMutsVsBranchesPairedPlotList) / 4), byrow=F)
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
#legend <- get_legend(numMutsInferredVsObsVsTruePairedPlotList[[1]] + theme(legend.box.margin=margin(0, 0, 0, 0), legend.position="bottom"))
#pGrid <- plot_grid(plotlist=lapply(numMutsInferredVsObsVsTruePairedPlotList, FUN=function(x) {x + theme(legend.position="none")}), align='vh', labels="AUTO", nrow=ceiling(length(numMutsInferredVsObsVsTruePairedPlotList) / 2), byrow=F)
pGrid <- plot_grid(plotlist=lapply(numMutsInferredVsObsVsTruePairedPlotList[1:4], FUN=function(x) {x + theme(legend.position="none")}), align='vh', labels="AUTO", ncol=1, byrow=T)
toSave <- plot_grid(pGrid, legend, ncol=1, rel_heights=c(1, 0.05))
outputFile <- paste0(gsub("/", "_", paramSets[1]), "-", gsub("/", "_", paramSets[length(paramSets)]), "_", filekeys[1], "-", filekeys[length(filekeys)], "_k", k, "_c", paste0(numCellsList, collapse="-c"), mutFilt, "_numMutsInferredVsObsVsTruePaired_A")
png(paste0(outputDir, "/", outputFile, ".png"), width=plotWidth*1.75, height=plotHeight*2, res=600, units="in"); plot(toSave); dev.off()
pGrid <- plot_grid(plotlist=lapply(numMutsInferredVsObsVsTruePairedPlotList[5:7], FUN=function(x) {x + theme(legend.position="none")}), align='vh', labels=c("E", "F", "G"), ncol=1, byrow=T)
toSave <- plot_grid(pGrid, legend, ncol=1, rel_heights=c(1, 0.05))
outputFile <- paste0(gsub("/", "_", paramSets[1]), "-", gsub("/", "_", paramSets[length(paramSets)]), "_", filekeys[1], "-", filekeys[length(filekeys)], "_k", k, "_c", paste0(numCellsList, collapse="-c"), mutFilt, "_numMutsInferredVsObsVsTruePaired_B")
png(paste0(outputDir, "/", outputFile, ".png"), width=plotWidth*1.75, height=plotHeight*2*3/4, res=600, units="in"); plot(toSave); dev.off()
#png(paste0(outputDir, "/", outputFile, ".png"), width=35, height=20, res=600, units="in"); plot(toSave); dev.off()

# pull out inferred vs obs muts for main figure
legend <- get_legend(numMutsInferredVsObsPairedPlotList[[1]] + theme(legend.box.margin=margin(0, 0, 0, 0), legend.position="bottom"))
pGrid <- plot_grid(plotlist=lapply(numMutsInferredVsObsPairedPlotList, FUN=function(x) {x + theme(legend.position="none")}), align='vh', labels="AUTO", ncol=2, byrow=T)
toSave <- plot_grid(pGrid, legend, ncol=1, rel_heights=c(1, 0.05))
outputFile <- paste0(gsub("/", "_", paramSets[1]), "-", gsub("/", "_", paramSets[length(paramSets)]), "_", filekeys[1], "-", filekeys[length(filekeys)], "_k", k, "_c", paste0(numCellsList, collapse="-c"), mutFilt, "_numMutsInferredVsObsPaired")
#png(paste0(outputDir, "/", outputFile, ".png"), width=35, height=20, res=600, units="in"); plot(toSave); dev.off()
png(paste0(outputDir, "/", outputFile, ".png"), width=plotWidth*1.5, height=plotHeight*2, res=600, units="in"); plot(toSave); dev.off()




###################################
## misc debugging
#v15 <- inferredPairedMutCountsList[[1]]
#v15_nan <- inferredPairedMutCountsList[[2]]
#
#v15_m <- melt(v15, id.vars=colnames(v15)[1:5])
#colnames(v15_m) <- c("left", "right", "treeBranch", "cell0", "cell1", "variable", "value")
#v15_c <- dcast(v15_m, left + right + cell0 + cell1 ~ variable + treeBranch)
#
#v15_nan_m <- melt(v15_nan, id.vars=colnames(v15_nan)[1:5])
#colnames(v15_nan_m) <- c("left", "right", "treeBranch", "cell0", "cell1", "variable", "value")
#v15_nan_c <- dcast(v15_nan_m, left + right + cell0 + cell1 ~ variable + treeBranch)
#
#ggplot(data=data.frame(v15=v15_c$inferred_sum, v15_nan=v15_nan_c$inferred_sum), aes(x=v15, y=v15_nan)) + geom_point() + geom_abline(slope=1, intercept=0, colour="red")
#ggplot(data=melt(data.frame(v15=v15_c$inferred_sum, v15_nan=v15_nan_c$inferred_sum)), aes(value, fill=variable)) + geom_histogram(position="dodge")
## ==> v15 is predicting too many muts compared to v15_nan; also confirmed with the numMutsInferredVsObsVsTruePaired plots
#
