library(devtools)
library(DESeq2)
load_all()
uploadDir <- system.file("extdata/bcbio", package = "bcbioRNASeq")
bcb <- loadRNASeq(
    uploadDir,
    interestingGroups = "group")
dds <- DESeqDataSetFromTximport(
    txi = txi(bcb),
    colData = colData(bcb),
    design = formula(~group)) %>%
    DESeq()
rld <- rlog(dds)
res <- results(dds)
use_data(bcb, dds, rld, res, compress = "xz", overwrite = TRUE)
