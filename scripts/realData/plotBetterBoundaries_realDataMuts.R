library(reshape2)
library(stringr)
library(ggplot2)
library(scales)
library(grid)
# Sun 24 Jul 2022 08:14:20 PM PDT
# script to find real data cells to plot with mutations

# exploration to find cell to plot: look at sseDat for high diff cells, then potentially look for regions of difference (run the stuff in calcIndvSconceSSE for cell of interest; then , plug in cell name, update arrows, plot
#colnames(sconceSconceMut) <- c("chr", "start", "end", "idx", "sconce", "sconceMut")
#colnames(sconce2SconceMut) <- c("chr", "start", "end", "idx", "sconce2", "sconceMut")
#sconceComps <- merge(sconceSconceMut, sconce2SconceMut, by=c("chr", "start", "end", "idx", "sconceMut"))
#sconceComps_sorted <- sconceComps[with(sconceComps, order(idx)),]
#sconceComps_sorted$sconceSconceMutDiff <- sconceComps_sorted$sconce - sconceComps_sorted$sconceMut
#sconceComps_sorted$sconce2SconceMutDiff <- sconceComps_sorted$sconce2 - sconceComps_sorted$sconceMut
#summary(subset(sconceComps_sorted, sconceSconceMutDiff !=0 & sconce2SconceMutDiff != 0))
#head(subset(sconceComps_sorted, abs(sconceSconceMutDiff) > 0.9 & abs(sconce2SconceMutDiff) > 0.5))

# find sse between (sconce, sconceMut), (sconce2 (mean), sconceMut (mean))
calcIndvSconceSSE <- function(sconceList, sconce2MeanFileList, sconceMutMeanFileList, cellIDs) {
  # cell | (sconce, sconceMut) | (sconce2, sconceMut)
  dat <- do.call(rbind, lapply(cellIDs, FUN=function(cellID) {
    sconceDat <- sconceList[[cellID]] #read.table(sconceFileList[grepl(cellID, sconceFileList, fixed=T)], stringsAsFactors=F, sep="\t", header=F)
    sconce2MeanDat <- read.table(sconce2MeanFileList[grepl(cellID, sconce2MeanFileList, fixed=T)], stringsAsFactors=F, sep="\t", header=F)
    sconceMutMeanDat <- read.table(sconceMutMeanFileList[grepl(cellID, sconceMutMeanFileList, fixed=T)], stringsAsFactors=F, sep="\t", header=F)

    colnames(sconce2MeanDat) <- colnames(sconceMutMeanDat)<- longHeader
    sconceDat$idx <- sconce2MeanDat$idx <- sconceMutMeanDat$idx <- 1:nrow(sconceDat)
  
    sconceSconceMut <- merge(sconceDat, sconceMutMeanDat, by=c("chr", "start", "end", "idx"))
    sconce2SconceMut <- merge(sconce2MeanDat, sconceMutMeanDat, by=c("chr", "start", "end", "idx"))
  
    sconceSconceMutSumSq <- sum((sconceSconceMut$copyNumber.x - sconceSconceMut$copyNumber.y)^2)
    sconce2SconceMutSumSq <- sum((sconce2SconceMut$copyNumber.x - sconce2SconceMut$copyNumber.y)^2)

    data.frame(cell=cellID, sconceSconceMut=sconceSconceMutSumSq, sconce2SconceMut=sconce2SconceMutSumSq)
  }))
  dat
}

calcPairsSconceSSE <- function(sconceList, pairsFileList, cellIDs, cellRegex) {
  # cells in pair | cell | variable (pair) | value (sumSq)
  dat <- do.call(rbind, lapply(pairsFileList, FUN=function(f) {
    splitFileName <- gsub(".bed", "", gsub("cov_unif_250kb", "", unlist(str_split(f, "__"))))
    pair <- gsub("pair_", "", splitFileName[2])
    cellID <- str_extract(splitFileName[3], cellRegex)
    pairDat <- read.table(f, stringsAsFactors=F, sep="\t", header=F)
    colnames(pairDat) <- longHeader
    pairDat$idx <- 1:nrow(pairDat)

    sconceDat <- sconceList[[cellID]]
    pairMerged <- merge(pairDat, sconceDat, by=c("chr", "start", "end", "idx"))
    sumSq <- sum((pairMerged$copyNumber.x - pairMerged$copyNumber.y)^2)
    data.frame(pair=pair, cell=cellID, variable="one pair", value=sumSq)
  }))
  dat
}

calcRealSSE <- function(paramSet, sconce2Key, sconceMutKey, sseFile, numCells, forceRecalc, cellRegex) {
  if(!forceRecalc && file.exists(paste0(sseFile, "_indv.txt"))) {
    indvDat <- read.table(paste0(sseFile, "_indv.txt"), header=T, sep="\t", stringsAsFactors=F)
    allPairDat <<- allPairDat
  } else {
    sconce2OutputFilePrefix <- paste0("output_", sconce2Key, "_", paramSet, "_k", k, "_c", numCells) # based on scAllP_*sh outBase variable
    sconceMutOutputFilePrefix <- paste0("output_", sconceMutKey, "_", paramSet, "_k", k, "*_c", numCells) # based on scAllP_*sh outBase variable # TODO change this key
    unfiltSconceFileList <- system(paste0("find ", dataDir, " -maxdepth 1 -name \"", sconceMutOutputFilePrefix, "*__sconce__*.bed\" | sort -V"), intern=T)
    unfiltSconce2MeanFileList <- system(paste0("find ", dataDir, " -maxdepth 1 -name \"", sconce2OutputFilePrefix, "*__mean.bed\" | sort -V"), intern=T)
    unfiltSconceMutMeanFileList <- system(paste0("find ", dataDir,     " -maxdepth 1 -name \"", sconceMutOutputFilePrefix, "*__mean.bed\" | sort -V"), intern=T)

    # if missing any type of output files to read, return null
    if(length(unfiltSconceFileList) == 0 || length(unfiltSconce2MeanFileList) == 0 || length(unfiltSconceMutMeanFileList) == 0) {
      return(NULL)
    }

    tumorDepthsFile <- paste0(dataDir, "/../tumor_cov_unif_250kb_", numCells)
    if(!file.exists(tumorDepthsFile)) {
      tumorDepthsFile <- paste0(dataDir, "/../tumor_cov_unif_250kb_short_", numCells)
    }
    tumorDepths <- read.table(tumorDepthsFile, stringsAsFactors=F)$V1
    cellIDs <- str_extract(tumorDepths, cellRegex)

    sconceFileList <- sapply(cellIDs, FUN=function(cellID) {unfiltSconceFileList[grepl(cellID, unfiltSconceFileList, fixed=T)]})
    sconce2MeanFileList <- sapply(cellIDs, FUN=function(cellID) {unfiltSconce2MeanFileList[grepl(cellID, unfiltSconce2MeanFileList, fixed=T)]})
    sconceMutMeanFileList <- sapply(cellIDs, FUN=function(cellID) {unfiltSconceMutMeanFileList[grepl(cellID, unfiltSconceMutMeanFileList, fixed=T)]})

    sconceList <<- lapply(sconceFileList, FUN=function(f) {
      currSconce <- read.table(f, stringsAsFactors=F, sep="\t", header=F)
      colnames(currSconce) <- longHeader
      currSconce$idx <- 1:nrow(currSconce)
      currSconce
    })

    indvDat <- calcIndvSconceSSE(sconceList, sconce2MeanFileList, sconceMutMeanFileList, cellIDs)
    write.table(indvDat, paste0(sseFile, "_indv.txt"), sep="\t", col.names=T, row.names=F, quote=F)
  }
  indvDat_m <- melt(indvDat)

  toPlot <- indvDat_m
  toPlot$variable <- factor(toPlot$variable, levels=c("sconceSconceMut", "sconce2SconceMut"))
  toPlot
}

changeSourceStripColors <- function(p, side, g=NULL) {
  if(is.null(g)) {
    g <- ggplot_gtable(ggplot_build(p))
  }
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

orderedSources <- c("sconce", "sconce2_mean", "sconceMut_mean")
cellTextCols <- c("gray10", "gray10")
sourceStripFills <- c("#9b59d0", SCONCE2="#e98686", SCONCEmut="#22a5e3")
sourceTextCols <- c("white", "gray10", "gray10")
names(sourceStripFills) <- orderedSources
longHeader <- c("chr", "start", "end", "copyNumber")

k <- 10
sconce2Key <- "scAllP_v26"
sconceMutKey <- "scAllPMuts_v16"

############
# navin data
numCellsList <- c(20)
dataDir <- "./Navin_Nature2011_hg19/sap_o_muts"
currOutputDir <- paste0(dataDir, "/plots/")
paramSet <- "Navin_Nature2011_hg19"
cellRegex <- "SRR[0-9]*."
if(!dir.exists(currOutputDir)) {
  dir.create(currOutputDir)
}
forceRecalc <- T

sseDatList <- list()
allPairDat <- NULL
for(numCells in numCellsList) {
  #for(key in filekeys) {
    datasetName <- paste0(paramSet, ", ", sconce2Key, ", ", sconceMutKey, ", ", numCells)
    print(paste0("reading SSE ", datasetName))
    sseFile <- paste0(currOutputDir, "output_", sconce2Key, "-", sconceMutKey, "_", paramSet, "_k", k, "_c", numCells, "_sumSqSconceOnePairMeanMedianMode") # based on scAllP_*sh outBase variable
    sseDat <- calcRealSSE(paramSet, sconce2Key, sconceMutKey, sseFile, numCells, forceRecalc, cellRegex)
    sseDatList[[datasetName]] <- sseDat

  #}
}

sseDat_m <- do.call(rbind, sseDatList)
sseDat <- dcast(sseDat_m, cell ~ variable, fun.aggregate=mean)

cellA <- "SRR053672."; minIdx <- 1540; maxIdx <- 1640; pointOfInterest <- 1584 # maybe ok! a couple spots where could draw arrows
# snps at this location, from ./dbsnpIntersection_tumor/SRR053672_withPooledMajorAllele.llr.dbsnp.minCells5_minAvgReads3.5.readsMat.llr_t3.bed
#    835 chr2    146591887       146591888       2       0
#    836 chr2    146591898       146591899       2       0
#    837 chr2    146591899       146591900       2       0
#    838 chr2    146591930       146591931       1       0
#    839 chr2    146615278       146615279       1       0

# copy numbers at this location
#        chr     start       end   idx sconceMut sconce  sconce2
# 5365   chr2 146250000 146500000  1584   3.00000      2 2.052600
# 5366   chr2 146500000 146750000  1585   3.00000      2 2.052600
# 5367   chr2 146750000 147000000  1586   3.00000      2 2.052600
# 5368   chr2 147000000 147250000  1587   3.00000      2 2.052600

# compare CN=2 vs CN=3 (number germline alleles)
# let omega = 10
# CN=2 ==>
# g=0, s=2; if N=2 ==> 0.0004772727   if N=1 ==> 0.005
# g=1, s=1;        ==> 0.2727273                 0.5
# g=2, s=0; f=0    ==> 0.9904773                 0.995

# CN=3 ==>
# g=0, s=3; if N=2 ==> 8.153409e-05   if N=1 ==> 0.005
# g=1, s=2;            0.05906278                0.335
# g=2, s=1;            0.333369                  0.665
# g=3, s=0;            0.9863503                 0.995

cellIDs <- c(cellA)
names(cellIDs) <- cellIDs
orderedCells <- cellIDs
cellID <- cellA

arrows <- data.frame(x=   c(rep(1586,3), rep(1598,3)),
                     xend=c(rep(1586,3), rep(1598,3)),
                     y=c(rep(5.5, 6)),
                     yend=c(rep(4.5,6)),
                     source=c("sconce", "sconce2_mean", "sconceMut_mean",       "sconce", "sconce2_mean", "sconceMut_mean"),
                     cell=c(rep(cellA,3), rep(cellA,3)))
arrows$source <- factor(arrows$source, levels=orderedSources)
arrows$cell <- factor(arrows$cell, levels=orderedCells)

sconce2OutputFilePrefix <- paste0("output_", sconce2Key, "_", paramSet, "_k", k)
sconceMutOutputFilePrefix <- paste0("output_", sconceMutKey, "_", paramSet, "_k", k)
unfiltSconceFileList <- system(paste0("find ", dataDir, " -maxdepth 1 -name \"", sconceMutOutputFilePrefix, "*__sconce__*.bed\" | sort -V"), intern=T)
unfiltSconce2MeanFileList <- system(paste0("find ", dataDir, " -maxdepth 1 -name \"", sconce2OutputFilePrefix, "*__mean.bed\" | sort -V"), intern=T)
unfiltSconceMutMeanFileList <- system(paste0("find ", dataDir,     " -maxdepth 1 -name \"", sconceMutOutputFilePrefix, "*__mean.bed\" | sort -V"), intern=T)

sconceFileList <- sapply(cellIDs, FUN=function(cellID) {unfiltSconceFileList[grepl(cellID, unfiltSconceFileList, fixed=T)]})
sconce2MeanFileList <- sapply(cellIDs, FUN=function(cellID) {unfiltSconce2MeanFileList[grepl(cellID, unfiltSconce2MeanFileList, fixed=T)]})
sconceMutMeanFileList <- sapply(cellIDs, FUN=function(cellID) {unfiltSconceMutMeanFileList[grepl(cellID, unfiltSconceMutMeanFileList, fixed=T)]})

readsDatA <- read.table(system(paste0("find ", dataDir, "/../ -name \"", cellA, "cov_unif_250kb\""), intern=T), stringsAsFactors=F, header=F)
sconceDatA <- read.table(sconceFileList[cellA])
sconce2MeanDatA <- read.table(sconce2MeanFileList[cellA])
sconceMutMeanDatA <- read.table(sconceMutMeanFileList[cellA])

dipAvg <- read.table(system(paste0("find ", dataDir, "/../diploid -name diploid_avg_cov_unif_250kb.bed"), intern=T), stringsAsFactors=F, header=F)

header <- longHeader[1:3]
colnames(readsDatA) <- c(header, "depth")
colnames(sconceDatA) <- c(header, "copyNumber")
colnames(sconce2MeanDatA) <- c(header, "copyNumber")
colnames(sconceMutMeanDatA) <- c(header, "copyNumber")
colnames(dipAvg) <- c(header, "mean", "var")

readsDatA$idx <- 1:nrow(readsDatA)
sconceDatA$idx <- 1:nrow(sconceDatA)
sconce2MeanDatA$idx <- 1:nrow(sconce2MeanDatA)
sconceMutMeanDatA$idx <- 1:nrow(sconce2MeanDatA)
dipAvg$idx <- 1:nrow(dipAvg)

sconceDatA$source <- "sconce"
sconce2MeanDatA$source <- "sconce2_mean"
sconceMutMeanDatA$source <- "sconceMut_mean"

readsDatA$cell <- sconceDatA$cell <- sconce2MeanDatA$cell <- sconceMutMeanDatA$cell <- factor(cellA, levels=orderedCells)
inferredPloidies <- rbind(sconceDatA, sconce2MeanDatA, sconceMutMeanDatA)

sconceDat <- sconceDatA
sconceDat$source <- NULL

facetLabels <- c("SCONCE", "SCONCE2 (mean)", "SCONCEmut (mean)", "cell A") # get all caps SCONCE and generic cell labels
names(facetLabels) <- c("sconce", "sconce2_mean", "sconceMut_mean", as.name(cellA))
inferredPloidies$source <- factor(inferredPloidies$source, levels=orderedSources)
inferredPloidies$cell <- factor(inferredPloidies$cell, levels=orderedCells)

scalingFactorsA <-  sapply(orderedSources, FUN=function(x) { sum(readsDatA$depth) / mean(inferredPloidies[inferredPloidies$source == x & inferredPloidies$cell == cellA, "copyNumber"], na.rm=T) / nrow(readsDatA)})
cnvScalingA <- mean(scalingFactorsA)
dipScaling <- mean(dipAvg$mean) / 2

copyNumberAxisMax <- ceiling(max(subset(readsDatA, idx > minIdx & idx < maxIdx)$depth)/cnvScalingA) + .5
p <- ggplot(subset(inferredPloidies, idx > minIdx & idx < maxIdx), aes(x=idx, y=copyNumber, colour=source)) +
  theme_bw() + guides(colour="none") +
  geom_line(data=subset(dipAvg, idx > minIdx & idx < maxIdx), aes(x=idx, y=mean / dipScaling), colour="lightblue", alpha=1) +
  geom_ribbon(data=subset(dipAvg, idx > minIdx & idx < maxIdx), aes(x=idx, ymin=mean/dipScaling - sqrt(var/(dipScaling^2)), ymax=mean/dipScaling + sqrt(var/(dipScaling^2))), fill="lightblue", alpha=0.3, inherit.aes=F) +
  geom_point(data=subset(readsDatA, idx > minIdx & idx < maxIdx), aes(x=idx, y=depth / cnvScalingA), alpha=0.5, size=0.85, colour="darkgray") +
  #geom_point(data=subset(readsDatA, idx > minIdx & idx < maxIdx), aes(x=idx, y=depth / cnvScalingA), alpha=1, size=2, colour="darkgray") +
  geom_step(size=1.15) +
  facet_wrap(~ source, labeller=as_labeller(facetLabels)) +
  geom_segment(data=arrows, aes(x=x, xend=xend, y=y, yend=yend), colour="black", arrow=arrow(length=unit(0.5, "char")), lineend="round", linejoin="round",  size=1.5) + # yellow arrows with black border
  geom_segment(data=arrows, aes(x=x, xend=xend, y=y, yend=yend), colour="#fff433", arrow=arrow(length=unit(0.5, "char")), lineend="round", linejoin="round",  size=.75) +
  scale_y_continuous(sec.axis=sec_axis(~ . * dipScaling, name="scaled read depth"), limits=c(0, copyNumberAxisMax), breaks=seq(0, copyNumberAxisMax, 2)) +
  labs(y="copy number", x="genomic index") +
  scale_colour_manual(values=sourceStripFills) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) 

plotWidth <- 8
plotHeight <- 4
plotRoot <- paste0(dataDir, "/plots/betterBoundsMuts_cells_", cellA, "_", paramSet, "_", sconce2Key, "-", sconceMutKey, "_k", k, "_c", numCells, "_", minIdx, "-", maxIdx)
outputFile <- paste0(plotRoot, "_all.png")
png(outputFile, width=plotWidth, height=plotHeight, res=600, units="in")
g <- changeSourceStripColors(p, "t")
grid::grid.draw(g)
dev.off()

outputFile <- paste0(plotRoot, "_all.pdf")
pdf(outputFile, width=plotWidth, height=plotHeight)
g <- changeSourceStripColors(p, "t")
grid::grid.draw(g)
dev.off()

