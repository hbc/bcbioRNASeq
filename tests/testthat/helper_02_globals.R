anno <- annotable(bcb)

sampleMetadataFile <- "sample_metadata.csv"
uploadDir <- system.file("extdata/bcbio", package = "bcbioRNASeq")

plotlist <- c("tree_row", "tree_col", "kmeans", "gtable")

lfc <- 0.25
