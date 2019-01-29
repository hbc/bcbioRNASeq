# sudo apt-get install libudunits2-dev
library(BiocManager)

BiocManager::install(c(
    "gtools",
    "hexbin",
    "ggplot2",
    "cowplot",
    "DESeq2",
    "DEGreport",
    "clusterProfiler",
    "org.Mm.eg.db",
    "hbc/bcbioRNASeq"
))

library(ggplot2)
library(cowplot)
library(DESeq2)
library(DEGreport)
library(bcbioRNASeq)

theme_set(
    theme_gray(base_size = 10)
)



load("data/bcb.rda")

lv = gtools::mixedsort(as.character(colData(bcb)[["sampleName"]]))
colData(bcb)[["sampleName"]] = factor(colData(bcb)[["sampleName"]], levels = lv)
lv = unique(gtools::mixedsort(as.character(colData(bcb)[["day"]])))
colData(bcb)[["day"]] = factor(colData(bcb)[["day"]], levels = lv)
interestingGroups(bcb) = "day"
saveData(bcb)

# Figure 1 ====
theme_update(
    legend.justification = "left",
    legend.position = "right"
)
pdf("figures/generalStats.pdf", width = 9, height = 7)
plot_grid(
    plotTotalReads(bcb),
    plotMappingRate(bcb),
    plotExonicMappingRate(bcb),
    plotIntronicMappingRate(bcb),
    plotGenesDetected(bcb),
    plotGeneSaturation(bcb),
    ncol = 2,
    labels = "AUTO",
    label_size = 12
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
pdf("figures/plotMeanSD.pdf", width = 8, height = 4)
plotMeanSD(bcb)
dev.off()



# Figure 4 ====
pdf("figures/plotDispEsts.pdf", width = 6, height = 6)
plotDispEsts(bcb)
dev.off()



# Figure 5 ====
pdf("figures/plotCorrelationHeatmap.pdf", width = 6, height = 6)
plotCorrelationHeatmap(bcb)
dev.off()



# Figure 6 ====
pdf("figures/plotPCA.pdf", width = 6, height = 4)
plotPCA(bcb, label = FALSE)
dev.off()



# Figure 7 ====
pdf("figures/plotPCACovariates.pdf", width = 8, height = 5)
plotPCACovariates(bcb)
dev.off()



dds <- as(bcb, "DESeqDataSet")
design(dds) <- ~day
dds <- DESeq(dds)

res <- results(
    dds,
    name = "day_7_vs_0",
    alpha = 0.05)


alphaSummary(dds, contrast = c("day", "7", "0"), alpha = c(0.1, 0.05))

# Figure 8 ====
pdf("figures/plotMA.pdf", width = 6, height = 6)
plotMeanAverage(res)
dev.off()



# Figure 9 ====
pdf("figures/plotVolcano.pdf", width = 6, height = 6)
plotVolcano(res)
dev.off()



# Figure 10 ====
pdf("figures/plotDEGHeatmap.pdf", width = 6, height = 6)
# DESeqResults, DESeqTransform
plotDEGHeatmap(res, counts = bcb, normalized = "vst")
dev.off()



# Figure 11 ====
pdf("figures/degPlot.pdf", width = 6, height = 3)
degPlot(
    object = bcb,
    res = res,
    n = 3,
    slot = "vst",
    log = FALSE,
    # Use the 'ann' argument for gene to symbol mapping
    ann = c("geneID", "geneName"),
    xs = "day")
dev.off()



dds_lrt <- DESeq(dds, test = "LRT", reduced = ~1)
res_lrt <- results(dds_lrt)
saveData(dds_lrt, res_lrt)



theme_set(
    theme_gray(base_size = 5)
)
ma <- counts(bcb, "rlog")[significants(res, fc = 2), ]
res_patterns <- degPatterns(
    ma,
    metadata = colData(bcb),
    time = "day",
    minc = 60)
saveData(dds_lrt, res_lrt, res_patterns)


# Figure 12 ====

pdf("figures/degPatterns.pdf", width = 8, height = 4)
res_patterns[["plot"]]
dev.off()



# Subset example
subset(res_patterns[["df"]], cluster == 8, select = "genes")

res_tbl <- resultsTables(res, lfc = 1)
topTables(res_tbl, n = 5)

# FA analysis
library(clusterProfiler)
sigGenes <- significants(res, fc = 1, padj = 0.05)
allGenes <- rownames(res)[!is.na(res[["padj"]])]
sigResults <- as.data.frame(res)[sigGenes, ]
foldChanges <- sigResults$log2FoldChange
names(foldChanges) <- rownames(sigResults)

ego <- enrichGO(
    gene = sigGenes,
    OrgDb = "org.Mm.eg.db",
    keyType = "ENSEMBL",
    ont = "BP",
    universe = allGenes,
    qvalueCutoff = 0.05,
    readable = TRUE
)

pdf("figures/dot_plot.pdf", width = 7, height = 5)
dotplot(ego, showCategory = 25)
dev.off()
pdf("figures/enrich_map.pdf", width = 6, height = 6)
emapplot(ego, n = 25)
dev.off()
pdf("figures/cnet_plot.pdf", width = 6, height = 6)
cnetplot(
    ego,
    categorySize = "pvalue",
    showCategory = 5,
    foldChange = foldChanges,
    node_label = FALSE,
)
dev.off()
