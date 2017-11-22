library(devtools)
library(DESeq2)
load_all()
bcb <- loadRNASeq(
    uploadDir = system.file("extdata/bcbio", package = "bcbioRNASeq"),
    interestingGroups = "group",
    ensemblVersion = 90)
dds <- DESeqDataSetFromTximport(
    txi = txi(bcb),
    colData = colData(bcb),
    design = formula(~group)) %>%
    DESeq()
rld <- rlog(dds)
res <- results(dds)
examples <- list(
    bcb = bcb,
    dds = dds,
    rld = rld,
    res = res
)
use_data(examples, compress = "xz", overwrite = TRUE)
