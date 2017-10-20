library(bcbioRNASeq)
library(DEGreport)
library(DESeq2)
library(ggplot2)
library(cowplot)

source("loadRNASeq.R")
loadData(bcb)

theme_set(
    theme_gray(base_size = 10)
)

# Figure 1 ====
theme_update(
    legend.justification = "left",
    legend.position = "right"
)
pdf("figures/generalStats.pdf", width = 6, height = 6)
plot_grid(
    plotTotalReads(bcb),
    plotMappingRate(bcb),
    plotExonicMappingRate(bcb),
    plotIntronicMappingRate(bcb),
    plotGenesDetected(bcb),
    plotGeneSaturation(bcb),
    ncol = 2,
    labels = "AUTO",
    label_size = 14
)
dev.off()



# Figure 2 ====
theme_update(
    legend.justification = "center",
    legend.position = "bottom"
)
pdf("figures/geneCounts.pdf", width = 8, height = 4)
plot_grid(
    plotCountsPerGene(bcb),
    plotCountDensity(bcb),
    ncol = 2,
    labels = "AUTO",
    label_size = 14
)
dev.off()



# Figure 3 ====
pdf("figures/plotMeanSD.pdf", width = 3, height = 9)
plotMeanSD(bcb)
dev.off()



# Figure 4 ====
pdf("figures/plotDispEsts.pdf", width = 6, height = 6)
plotDispEsts(bcb)
dev.off()



# Figure 5 ====
pdf("figures/plotCorrelationHeatmap.pdf", width = 6, height = 4)
plotCorrelationHeatmap(bcb)
dev.off()



# Figure 6 ====
pdf("figures/plotPCA.pdf", width = 6, height = 4)
plotPCA(bcb, label = FALSE)
dev.off()



# Figure 7 ====
pdf("figures/plotPCACovariates.pdf", width = 6, height = 4.5)
plotPCACovariates(bcb)
dev.off()



dds <- DESeqDataSetFromTximport(
     txi = bcbio(bcb, "tximport"),
     colData = colData(bcb),
     design = formula(~group)) %>%
     DESeq()
res <- results(
    dds,
    contrast = c(
        factor = "group",
        numerator = "day7",
        denominator = "normal"),
    alpha = 0.05)
saveData(dds, res)



# Figure 8 ====
pdf("figures/plotMA.pdf", width = 6, height = 6)
plotMA(res)
dev.off()



# Figure 9 ====
pdf("figures/plotVolcano.pdf", width = 6, height = 6)
plotVolcano(res)
dev.off()



# Figure 10 ====
pdf("figures/plotDEGHeatmap.pdf", width = 6, height = 6)
# DESeqResults, DESeqTransform
plotDEGHeatmap(res, counts = assays(bcb)[["rlog"]])
dev.off()



# Figure 11 ====
pdf("figures/degPlot.pdf", width = 6, height = 4)
degPlot(
    bcb,
    res = res,
    n = 3,
    slot = "normalized",
    # Use the 'ann' argument for gene to symbol mapping
    ann = c("ensgene", "symbol"),
    xs = "group")
dev.off()



ddsLRT <- DESeqDataSetFromTximport(
    txi = bcbio(bcb, "tximport"),
    colData = colData(bcb),
    design = formula(~group)) %>%
    DESeq()
resLRT <- results(ddsLRT)
saveData(ddsLRT, resLRT)



# Use `consensusCluster = FALSE` for faster initial visualization
# For `consensusCluster = TRUE`, run on O2
resPatterns <- degPatterns(
    counts(bcb, "rlog")[significants(res), ],
    metadata = colData(bcb),
    time = "group",
    col = NULL,
    concensusCluster = TRUE)
saveData(resPatterns)



# Subset example
subset(resPatterns[["df"]], cluster == 8, select = "genes")



# Figure 12 ====
pdf("figures/degPatterns.pdf", width = 6, height = 4)
resPatterns[["plot"]]
dev.off()
