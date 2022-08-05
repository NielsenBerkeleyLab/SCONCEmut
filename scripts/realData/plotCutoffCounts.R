library(ggplot2)
library(reshape2)
library(cowplot)

# usage:
# Rscript plotCutoffCounts.R dbsnpIntersection/pooledDiploidMajorAllele.llr.dbsnp.countLlrSummary.cutoffCounts.txt
args <- commandArgs(trailingOnly=TRUE)

tab <- read.table(args[1], sep="\t", header=F)
colnames(tab) <- c("llrCutoff", "numInDbsnp", "numNotInDbsnp", "fracInDbsnp", "fracNotInDbsnp")
tab <- subset(tab, numInDbsnp + numNotInDbsnp > 0)

tab_num <- tab[,c(1,2,3)]
tab_frac <- tab[,c(1,4,5)]

tab_num_m <- melt(tab_num, id.vars="llrCutoff")
tab_frac_m <- melt(tab_frac, id.vars="llrCutoff")

tab_num_m$variable <- factor(tab_num_m$variable, levels=c("numNotInDbsnp", "numInDbsnp"))
tab_frac_m$variable <- factor(tab_frac_m$variable, levels=c("fracNotInDbsnp", "fracInDbsnp"))

cutoff.95 <- tab_frac[min(which(tab_frac$fracInDbsnp >= 0.95)), "llrCutoff"]
cutoff1 <- tab_frac[min(which(tab_frac$fracInDbsnp >= 1)), "llrCutoff"]

fillColors <- c("#e98686", "#22a5e3") # pink, blue

p <- ggplot(tab_frac_m, aes(x=llrCutoff, y=value, fill=variable)) + geom_bar(stat="identity", position="stack") + geom_text(data=tab_num_m, aes(x=llrCutoff, y=ifelse(variable == "numNotInDbsnp", 1, 0), label=value, angle=90, hjust=ifelse(variable == "numNotInDbsnp", 1, 0)), colour="black", inherit.aes=F) + theme_bw() + labs(x="loglikelihood ratio cutoff", y="fraction of diploid sites in dbSNP above LLR cutoff", fill="dbSNP\nmembership") + scale_fill_manual(labels=c("not in dbSNP", "in dbSNP"), values=fillColors) + geom_vline(xintercept=cutoff.95, color="blue")
ggsave(paste0(args[1], "_llrCutoffCounts_labeled.png"), width=6, height=5)

fracPlot <- ggplot(tab_frac_m, aes(x=llrCutoff, y=value, fill=variable)) + geom_bar(stat="identity", position="stack") + theme_bw() + labs(x="loglikelihood ratio cutoff", y="fraction of diploid sites in dbSNP", fill="dbSNP\nmembership") + scale_fill_manual(labels=c("not in dbSNP", "in dbSNP"), values=fillColors) + geom_vline(xintercept=cutoff.95, color="blue")
ggsave(paste0(args[1], "_llrCutoffCounts.png"), width=6, height=5)

print(paste0("95% llr cutoff is: ",cutoff.95))

fracPlot <- ggplot(subset(tab_frac_m, llrCutoff <= cutoff1), aes(x=llrCutoff, y=value, fill=variable)) + geom_bar(stat="identity", position="stack") + theme_bw() + labs(x="loglikelihood ratio cutoff", y="fraction of sites in dbSNP\nabove LLR cutoff", fill="dbSNP membership") + scale_fill_manual(labels=c("not in dbSNP", "in dbSNP"), values=fillColors) + geom_vline(xintercept=cutoff.95, color="blue")

countPlot <- ggplot(subset(tab_num_m, llrCutoff >= 1 & llrCutoff <=cutoff1), aes(x=llrCutoff, y=value, fill=variable)) + geom_bar(stat="identity", position="stack") + theme_bw() + labs(x="loglikelihood ratio cutoff", y="number of sites\nabove LLR cutoff") + geom_vline(xintercept=cutoff.95, colour="blue") + scale_y_continuous(labels = scales::scientific) + scale_fill_manual(labels=c("not in dbSNP", "in dbSNP"), values=fillColors)

plotList <- list(fracPlot, countPlot)
legend <- get_legend(plotList[[1]] + theme(legend.box.margin=margin(0, 0, 0, 0), legend.position="bottom") + guides(colour=guide_legend(nrow=1)))
pGrid <- plot_grid(plotlist=lapply(plotList, FUN=function(x) {x + theme(legend.position="none")}), align='vh', labels="AUTO", nrow=2, byrow=F)
toSave <- plot_grid(pGrid, legend, ncol=1, rel_heights=c(1, 0.05))
 
png(paste0(args[1], "_llrCutoffCounts_fracCount.png"), width=8, height=5.5, res=600, units="in"); plot(toSave); dev.off()
pdf(paste0(args[1], "_llrCutoffCounts_fracCount.pdf"), width=8, height=5.5); plot(toSave); dev.off()

