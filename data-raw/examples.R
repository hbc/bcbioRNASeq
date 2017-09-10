library(devtools)
library(DESeq2)
load_all()
extraDir <- system.file("extra", package = "bcbioRNASeq")
uploadDir <- file.path(extraDir, "bcbio")
bcb <- loadRNASeqRun(uploadDir)
dds <- DESeqDataSetFromTximport(
    txi = txi(bcb),
    colData = colData(bcb),
    design = formula(~group)) %>%
    DESeq
rld <- rlog(dds)
res <- results(dds)
use_data(bcb, dds, rld, res, compress = "xz", overwrite = TRUE)
