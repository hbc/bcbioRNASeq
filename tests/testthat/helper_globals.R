# nolint start

# bcbio default color palette support
# options(
#     bcbio.discrete.color = ggplot2::scale_color_viridis_d(),
#     bcbio.discrete.fill = ggplot2::scale_fill_viridis_d()
# )

uploadDir <- system.file("extdata/bcbio", package = "bcbioRNASeq")
lfc <- 0.25

dds_small <- as(deseq_small, "DESeqDataSet")
vst_small <- as(deseq_small, "DESeqTransform")
# Consider using a pre-saved object here, for speed.
rlog_small <- DESeq2::rlog(dds_small)
res_small <- as(deseq_small, "DESeqResults")
res_tables <- DESeqResultsTables(res_small)

# nolint end
