# nolint start

# bcbio default color palette support
options(
    bcbio.discrete.color = ggplot2::scale_color_viridis_d(),
    bcbio.discrete.fill = ggplot2::scale_fill_viridis_d()
)

lfc <- 0.25
pheatmapNames <- c("tree_row", "tree_col", "kmeans", "gtable")

# DESeqTransform
rld_small <- DESeq2::rlog(dds_small)
vst_small <- DESeq2::varianceStabilizingTransformation(dds_small)

# nolint end
