load("../kidney.rda")

library(cowplot)

library(bcbioRNASeq)
library(DESeq2)
library(ggplot2)
library(DEGreport)
design <- formula(~1)

basefontsize=10

theme_set(theme_light(base_size=basefontsize))

pdf("figure_generalstats.pdf", height=9, width=7.5)
# Figure 1
plot_grid(
  plotTotalReads(bcb),
  plotMappingRate(bcb),
  plotExonicMappingRate(bcb),
  plotIntronicMappingRate(bcb),
  plotGenesDetected(bcb),
  plotGeneSaturation(bcb),
  nrow=3,
  labels=c("A", "B", "C", "D", "E", "F"))
dev.off()

pdf("figure_vst.pdf", height=7.5, width=5)
# Figure 2
plotMeanSD(bcb)
dev.off()

pdf("figure_dispersions.pdf", height=4.5, width=7.5)
#figure3
plotDispEsts(bcb)
dev.off()

pdf("figure_genecounts.pdf", height=3, width=7.5)
plot_grid(
  plotCountsPerGene(bcb),
  plotCountDensity(bcb),
  nrow=1,
  labels=c("A", "B"))
dev.off()

pdf("figure_pca.pdf", height=3, width=7.5)
plot_grid(
  plotPCA(bcb, label=FALSE),
  plotPCA(bcb, label=TRUE),
  nrow=1,
  labels=c("A", "B"))
dev.off()

pdf("figure_correlationheatmap.pdf", height=4.5, width=6.5)
plotCorrelationHeatmap(bcb)
dev.off()

pdf("figure_pca_covariates.pdf", height=6.5, width=6.5)
cols2show = metrics(bcb) %>% colnames() %>% setdiff(., "x53Bias")
plotPCACovariates(bcb, metrics = cols2show, fdr=0.2)
dev.off()


design <- formula(~group)
alpha <- 0.05
lfc <- 1
txobject <- bcbio(bcb, "tximport")
colData <- colData(bcb)
colData(bcb)[["group"]] <- relevel(colData(bcb)[["group"]], "normal")

ddsde <- DESeqDataSetFromTximport(txi=txobject, colData=colData(bcb), design=design) %>% DESeq()

alphaSummary(ddsde, name = "group_day7_vs_normal")
res <- results(ddsde, name = "group_day7_vs_normal", alpha=0.05)

pdf("figure_ma.pdf", height=4.5, width=7.5)
plotMA(res)
dev.off()


pdf("figure_volcano.pdf", height=7.5, width=7.5)
plotVolcano(res, shadeColor = "green", pointColor = "grey66",
            pointOutlineColor = "grey44", pointAlpha = 0.8, shadeAlpha = 0.20)
dev.off()


pdf("figure_degheatmap.pdf", height=4.5, width=7.5)
plotDEGHeatmap(res, assays(bcb)[["rlog"]])
dev.off()

pdf("figure_deggenes.pdf", height=4.5, width=7.5)
degPlot(bcb, res = res, n = 3, xs = "group", slot = "normalized", ann = c("ensgene", "symbol"))
dev.off()

ddsdetime <- DESeqDataSetFromTximport(txi=txobject, colData=colData(bcb), design=design) %>% DESeq(test="LRT", reduced=~1)

sign <- row.names(res)[order(res$padj)[1:500]]
res_pattern <- degPatterns(counts(bcb, "rlog")[sign,],colData(bcb), time="group", col=NULL)

pdf("figure_patterns.pdf", height=4.5, width=7.5)
res_pattern[["plot"]]
dev.off()
