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
#numCellsList <- c(20, 40)
#numCellsList <- c(20, 40, 60, 80, 100, 120)
filekeys <- c(
              #"scAllP_v21")
              #"scAllP_v24")
              #"scAllP_v24_nearest10")
              "scAllPMuts_v16", "scAllP_v26")
#paramSets <- c("muts/params14", "muts/params15", "muts/params16", "muts/params17")
paramSets <- c("muts/params14", "muts/params15", "muts/params16", "muts/params17", "muts/params32", "muts/params33", "muts/params34") #, "muts/params20")
#paramSets <- c("muts/params33", "muts/params34") #, "muts/params20")
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

#print("plotting giant rf dist plot")
#legend <- get_legend(rfPlotList[[1]] + theme(legend.box.margin=margin(0, 0, 0, 0), legend.position="bottom") + guides(colour=guide_legend(nrow=1)))
#pGrid <- plot_grid(plotlist=lapply(rfPlotList, FUN=function(x) {x + theme(legend.position="none")}), align='vh', labels="AUTO", nrow=4, byrow=F)
#toSave <- plot_grid(pGrid, legend, ncol=1, rel_heights=c(1, 0.05))
#outputFile <- paste0(gsub("/", "_", paramSets[1]), "-", gsub("/", "_", paramSets[length(paramSets)]), "_", filekeys[1], "-", filekeys[length(filekeys)], "_k", k, "_c", paste0(numCellsList, collapse="-c"), "_rfDistWithZZS")
#png(paste0(outputDir, "/", outputFile, ".png"), width=25, height=20, res=600, units="in")
#plot(toSave); dev.off()


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
#pdf(paste0(outputDir, "/", outputFile, ".pdf"), width=plotWidth, height=plotHeight); plot(toSave); dev.off()





## faceted plot: sections/facets are distance metric, true/sconce/mean/median/mode[t2+t3] are box plots
#print("plotting faceted plot")
#for(paramSet in paramSets) {
#  p <- makeFacetedRFdistPlot(combinedTreeDistList[[paramSet]], NA)
#  p <- p + theme(strip.text.x = element_text(size = 6)) # shrink facet labels
#  combinedTreeDistFacetPlotList[[paramSet]] <- p
#}
#legend <- get_legend(combinedTreeDistFacetPlotList[[1]] + theme(legend.box.margin=margin(0, 0, 0, 0), legend.position="bottom") + guides(colour=guide_legend(nrow=1)))
#pGrid <- plot_grid(plotlist=lapply(combinedTreeDistFacetPlotList, FUN=function(x) {x + theme(legend.position="none")}), align='vh', labels="AUTO", nrow=2)
#toSave <- plot_grid(pGrid, legend, ncol=1, rel_heights=c(1, 0.05))
#outputFile <- paste0(gsub("/", "_", paramSets[1]), "-", gsub("/", "_", paramSets[length(paramSets)]), "_", filekeys[1], "-", filekeys[length(filekeys)], "_k", k, "_c", paste0(numCellsList, collapse="-c"), "_combinedRfDistWithZZS_facet")
##png(paste0(outputDir, "/", outputFile, ".png"), width=15, height=10, res=600, units="in")
##plot(toSave); dev.off()
#
#png(paste0(outputDir, "/", outputFile, ".png"), width=plotWidth, height=plotHeight, res=600, units="in"); plot(toSave); dev.off()
#pdf(paste0(outputDir, "/", outputFile, ".pdf"), width=plotWidth, height=plotHeight); plot(toSave); dev.off()
##save_plot(paste0(outputDir, "/", outputFile, ".eps"), toSave, device=cairo_ps, dpi=600, base_width=plotWidth, base_height=plotHeight)
#
#
#
#
#
#
## plot combining each subset of cells within a paramset
#if(length(numCellsList) > 1) {
#  print("combining cell subsets across paramsets")
#  combinedTreeDistList <- list()
#  combinedTreeDistPlotList <- list()
#  combinedTreeDistFiltPlotList <- list()
#  combinedTreeDistFacetPlotList <- list()
#  print("plotting cell subsets across paramsets")
#  for(paramSet in paramSets) {
#    combinedTreeDistList[[paramSet]] <- do.call(rbind, treeDistList[grepl(gsub("muts/params", "mp", paramSet), names(treeDistList))])
#    combinedTreeDistPlotList[[paramSet]] <- makeRFdistPlot(combinedTreeDistList[[paramSet]], NA)
#  }
#  legend <- get_legend(combinedTreeDistPlotList[[1]] + theme(legend.box.margin=margin(0, 0, 0, 0), legend.position="bottom") + guides(colour=guide_legend(nrow=1)))
#  pGrid <- plot_grid(plotlist=lapply(combinedTreeDistPlotList, FUN=function(x) {x + theme(legend.position="none")}), align='vh', labels="AUTO", nrow=2)
#  toSave <- plot_grid(pGrid, legend, ncol=1, rel_heights=c(1, 0.05))
#  outputFile <- paste0(gsub("/", "_", paramSets[1]), "-", gsub("/", "_", paramSets[length(paramSets)]), "_", filekeys[1], "-", filekeys[length(filekeys)], "_k", k, "_c", paste0(numCellsList, collapse="-c"), "_combinedRfDistWithZZS")
#  #png(paste0(outputDir, "/", outputFile, ".png"), width=15, height=10, res=600, units="in")
#  #plot(toSave); dev.off()
#
#  #png(paste0(outputDir, "/", outputFile, ".png"), width=plotWidth, height=plotHeight, res=600, units="in"); plot(toSave); dev.off()
#  #save_plot(paste0(outputDir, "/", outputFile, ".eps"), toSave, device=cairo_ps, dpi=600, base_width=plotWidth, base_height=plotHeight)
#
#  # small plot
#  # euc, cnp2, zzs on sconce, t2+t3
#  print("plotting just euc/cnp/zzs on sconce, t2+t3")
#  for(paramSet in paramSets) {
#    combinedTreeDistFiltPlotList[[paramSet]] <- makeRFdistPlot(subset(combinedTreeDistList[[paramSet]], variable %in% c("eucSconce", "cnpSconce", "zzsSconce", "t2_t3")), NA)
#  }
#  legend <- get_legend(combinedTreeDistFiltPlotList[[1]] + theme(legend.box.margin=margin(0, 0, 0, 0), legend.position="bottom") + guides(colour=guide_legend(nrow=1)))
#  pGrid <- plot_grid(plotlist=lapply(combinedTreeDistFiltPlotList, FUN=function(x) {x + theme(legend.position="none")}), align='vh', labels="AUTO", nrow=2)
#  toSave <- plot_grid(pGrid, legend, ncol=1, rel_heights=c(1, 0.05))
#  outputFile <- paste0(gsub("/", "_", paramSets[1]), "-", gsub("/", "_", paramSets[length(paramSets)]), "_", filekeys[1], "-", filekeys[length(filekeys)], "_k", k, "_c", paste0(numCellsList, collapse="-c"), "_combinedRfDistWithZZS_filt")
#  #png(paste0(outputDir, "/", outputFile, ".png"), width=15, height=10, res=600, units="in")
#  #plot(toSave); dev.off()
#
#  png(paste0(outputDir, "/", outputFile, ".png"), width=plotWidth, height=plotHeight, res=600, units="in"); plot(toSave); dev.off()
#  pdf(paste0(outputDir, "/", outputFile, ".pdf"), width=plotWidth, height=plotHeight); plot(toSave); dev.off()
#  #save_plot(paste0(outputDir, "/", outputFile, ".eps"), toSave, device=cairo_ps, dpi=600, base_width=plotWidth, base_height=plotHeight)
#
#
#  # faceted plot: sections/facets are distance metric, true/sconce/mean/median/mode[t2+t3] are box plots
#  print("plotting faceted plot")
#  for(paramSet in paramSets) {
#    p <- makeFacetedRFdistPlot(combinedTreeDistList[[paramSet]], NA)
#    p <- p + theme(strip.text.x = element_text(size = 6)) # shrink facet labels
#    combinedTreeDistFacetPlotList[[paramSet]] <- p
#  }
#  legend <- get_legend(combinedTreeDistFacetPlotList[[1]] + theme(legend.box.margin=margin(0, 0, 0, 0), legend.position="bottom") + guides(colour=guide_legend(nrow=1)))
#  pGrid <- plot_grid(plotlist=lapply(combinedTreeDistFacetPlotList, FUN=function(x) {x + theme(legend.position="none")}), align='vh', labels="AUTO", nrow=2)
#  toSave <- plot_grid(pGrid, legend, ncol=1, rel_heights=c(1, 0.05))
#  outputFile <- paste0(gsub("/", "_", paramSets[1]), "-", gsub("/", "_", paramSets[length(paramSets)]), "_", filekeys[1], "-", filekeys[length(filekeys)], "_k", k, "_c", paste0(numCellsList, collapse="-c"), "_combinedRfDistWithZZS_facet")
#  #png(paste0(outputDir, "/", outputFile, ".png"), width=15, height=10, res=600, units="in")
#  #plot(toSave); dev.off()
#
#  png(paste0(outputDir, "/", outputFile, ".png"), width=plotWidth, height=plotHeight, res=600, units="in"); plot(toSave); dev.off()
#  pdf(paste0(outputDir, "/", outputFile, ".pdf"), width=plotWidth, height=plotHeight); plot(toSave); dev.off()
#  #save_plot(paste0(outputDir, "/", outputFile, ".eps"), toSave, device=cairo_ps, dpi=600, base_width=plotWidth, base_height=plotHeight)
#
#
#  #allDists <- do.call(rbind, combinedTreeDistList)
#  #outputFile <- paste0(gsub("/", "_", paramSets[1]), "-", gsub("/", "_", paramSets[length(paramSets)]), "_", filekeys[1], "-", filekeys[length(filekeys)], "_k", k, "_c", paste0(numCellsList, collapse="-c"), "_combinedRfDistWithZZS")
#  #writeMedianRFdistTexFiles(allDists, paste0(outputDir, "/", outputFile))
#}
#
#
