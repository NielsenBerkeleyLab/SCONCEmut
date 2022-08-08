# Tue 31 May 2022 03:26:26 PM PDT
# script to plot SSE and breakpoint detection accuracy on paired data, with mutations. Based heavily on plotSSEandBreakpointPairs.R

# will create countedBreakpointsTable*txt, *sumSqSconceOnePairMeanMedianMode_indv.txt, and *sumSqSconceOnePairMeanMedianMode_pairs.txt if they don't exist, and will read from file if they do

# arg 1: k
# arg 2: numCellsList (ex 20,40)
# arg 3: filekeys (ex scAllP_v21,scAllP_v23)
# arg 4+: paramSets (ex "muts/params21" "muts/params22")

source("/space/s1/sandra/src/input/treeSim/scripts/readBedFilesPairs.R")

k <- 10
numCellsList <- 20
#numCellsList <- c(20, 40, 60, 80, 100, 120)
filekeys <- c(
              "scAllPMuts_v16", "scAllP_v26")
paramSets <- c("muts/params14", "muts/params15", "muts/params16", "muts/params17", "muts/params32", "muts/params33", "muts/params34")
forceRecalc <- T
inclAneu <- F

args <- commandArgs(trailingOnly=TRUE)
if(length(args) > 0) {
  k <- args[1]
  numCellsList <- unlist(strsplit(args[2], ",")) # 20,40,60,80,100,120
  filekeys <- unlist(strsplit(args[3], ",")) # scAllP_v21
  paramSets <- args[4:length(args)] # "muts/params21" "muts/params22"
}

breakpointDatList <- list()
breakpointPlotList <- list()
sseDatList <- list()
ssePlotList <- list()

for(numCells in numCellsList) {
  for(paramSet in paramSets) {
    for(key in filekeys) {
      mutFilt <- ""
      currOutputDir <- paste0(dataDir, "/", paramSet, "/plots/")
      if(!dir.exists(currOutputDir)) {
        dir.create(currOutputDir)
      }
      datasetName <- paste0(paramSet, ", ", key, ", ", numCells, mutFilt)
      print(paste0("reading SSE ", datasetName))
      sseFile <- paste0(currOutputDir, "output_", key, "_", gsub("/", "_", paramSet), "_k", k, "_c", numCells, mutFilt, "_sumSqSconceOnePairMeanMedianMode") # based on scAllP_*sh outBase variable
      sseDat <- calcAllSSE(paramSet, key, sseFile, numCells, mutFilt=mutFilt, forceRecalc=forceRecalc, inclAneu=inclAneu)
      if(is.null(sseDat)) {
        next
      }
      sseDatList[[datasetName]] <- sseDat

      print(paste0("reading breakpoint ", datasetName))
      breakpointFile <- paste0(currOutputDir, "output_", key, "_", gsub("/", "_", paramSet), "_k", k, "_c", numCells, mutFilt, "_breakpointDistSconceOnePairMeanMedianMode") # based on scAllP_*sh outBase variable
      breakpointDat <- calcAllBreakpoints(paramSet, key, breakpointFile, numCells, mutFilt=mutFilt, forceRecalc=forceRecalc, inclAneu=inclAneu)
      breakpointDatList[[datasetName]] <- breakpointDat

      print(paste0("saving tex files ", datasetName))
      medianDistOmegaTexFile <- paste0(currOutputDir, "output_", key, "_", gsub("/", "_", paramSet), "_k", k, "_c", numCells, mutFilt, "_medianDistOmegaSconceOnePairMeanMedianMode") # based on scAllP_*sh outBase variable
      breakpointDat$paramSet <- paramSet
      writeMedianDistOmegaTexFiles(breakpointDat, medianDistOmegaTexFile)
    }
  }
  for(paramSet in paramSets) {
    for(key in filekeys) {
      datasetName <- paste0(paramSet, ", ", key, ", ", numCells, mutFilt)
      shortName <- paste0(gsub("muts/params", "mp", paramSet), "_", key, "_c", numCells, mutFilt)
      if(datasetName %in% names(sseDatList)) {
        print(paste0("plotting SSE ", datasetName))
        ssePlot <- makeSSEPlot(sseDatList[[datasetName]], shortName)
        ssePlotList[[datasetName]] <- ssePlot
      }
      if(datasetName %in% names(breakpointDatList)) {
        print(paste0("plotting breakpoint ", datasetName))
        breakpointPlot <- makeBreakpointPlot(breakpointDatList[[datasetName]], shortName)
        breakpointPlotList[[datasetName]] <- breakpointPlot
      }
    }
  }
}

# make sconce2 vs sconceMut sse plots
sconce2MutSSEDatList <- list()
sconce2MutSSEPlotList <- list()
for(paramSet in paramSets) {
  sseDat <- sseDatList[grepl(paramSet, names(sseDatList))]
  sseDat[grepl("Mut", names(sseDat))][[1]]$program <- "SCONCEmut"
  sseDat[!grepl("Mut", names(sseDat))][[1]]$program <- "SCONCE2"
  sseDat <- rbind(sseDat[[1]], sseDat[[2]])
  sseDat$program <- factor(sseDat$program)
  sconce2MutSSEDatList[[paramSet]] <- sseDat
  sconce2MutSSEPlotList[[paramSet]] <- makeSSEPlotCompareMuts(sseDat, NA)
}
legend <- get_legend(sconce2MutSSEPlotList[[1]] + theme(legend.box.margin=margin(0, 0, 0, 0), legend.position="bottom", legend.title=element_blank()))
pGrid <- plot_grid(plotlist=lapply(sconce2MutSSEPlotList, FUN=function(x) {x + theme(legend.position="none")}), align='vh', labels="AUTO", ncol=2, byrow=T)
toSave <- plot_grid(pGrid, legend, ncol=1, rel_heights=c(1, 0.05))
outputFile <- paste0(gsub("/", "_", paramSets[1]), "-", gsub("/", "_", paramSets[length(paramSets)]), "_", filekeys[1], "-", filekeys[length(filekeys)], "_k", k, "_c", paste0(numCellsList, collapse="-c"), "_SSESconce2SconceMut")
png(paste0(outputDir, "/", outputFile, ".png"), width=plotWidth, height=plotHeight*1.5, res=600, units="in")
plot(toSave); dev.off()

# make sconce2 vs sconceMut breakpoint plots
sconce2MutBreakpointDatList <- list()
sconce2MutBreakpointPlotList <- list()
for(paramSet in paramSets) {
  breakpointDat <- breakpointDatList[grepl(paramSet, names(breakpointDatList))]
  breakpointDat[grepl("Mut", names(breakpointDat))][[1]]$program <- "SCONCEmut"
  breakpointDat[!grepl("Mut", names(breakpointDat))][[1]]$program <- "SCONCE2"
  breakpointDat <- rbind(breakpointDat[[1]], breakpointDat[[2]])
  breakpointDat$program <- factor(breakpointDat$program)
  sconce2MutBreakpointDatList[[paramSet]] <- breakpointDat
  sconce2MutBreakpointPlotList[[paramSet]] <- makeBreakpointPlotCompareMuts(breakpointDat, NA)
}
legend <- get_legend(sconce2MutBreakpointPlotList[[1]] + theme(legend.box.margin=margin(0, 0, 0, 0), legend.position="bottom", legend.title=element_blank()))
pGrid <- plot_grid(plotlist=lapply(sconce2MutBreakpointPlotList, FUN=function(x) {x + theme(legend.position="none")}), align='vh', labels="AUTO", ncol=2, byrow=T)
toSave <- plot_grid(pGrid, legend, ncol=1, rel_heights=c(1, 0.05))
outputFile <- paste0(gsub("/", "_", paramSets[1]), "-", gsub("/", "_", paramSets[length(paramSets)]), "_", filekeys[1], "-", filekeys[length(filekeys)], "_k", k, "_c", paste0(numCellsList, collapse="-c"), "_breakpointSconce2SconceMut")
png(paste0(outputDir, "/", outputFile, ".png"), width=plotWidth, height=plotHeight*1.6, res=600, units="in")
plot(toSave); dev.off()

allBreakpointDat <- do.call(rbind, breakpointDatList[grepl("Mut", names(breakpointDatList))])
allBreakpointDat$paramSet <- sapply(rownames(allBreakpointDat), FUN=function(x) {unlist(str_split(x, ","))[1]})
outputFile <- paste0(outputDir, "/", gsub("/", "_", paramSets[1]), "-", gsub("/", "_", paramSets[length(paramSets)]), "_", filekeys[1], "-", filekeys[length(filekeys)], "_k", k, "_c", paste0(numCellsList, collapse="-c"), mutFilt, "_combinedMedianDistOmegaSconceMut")
writeMedianDistOmegaTexFiles(allBreakpointDat, outputFile)

allBreakpointDat <- do.call(rbind, breakpointDatList[!grepl("Mut", names(breakpointDatList))])
allBreakpointDat$paramSet <- sapply(rownames(allBreakpointDat), FUN=function(x) {unlist(str_split(x, ","))[1]})
outputFile <- paste0(outputDir, "/", gsub("/", "_", paramSets[1]), "-", gsub("/", "_", paramSets[length(paramSets)]), "_", filekeys[1], "-", filekeys[length(filekeys)], "_k", k, "_c", paste0(numCellsList, collapse="-c"), mutFilt, "_combinedMedianDistOmegaSconce2")
writeMedianDistOmegaTexFiles(allBreakpointDat, outputFile)

