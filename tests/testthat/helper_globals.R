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
res_small <- as(deseq_small, "DESeqResults")

# nolint end
