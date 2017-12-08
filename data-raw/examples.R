devtools::load_all()

extdataDir <- file.path("inst", "extdata")
uploadDir <- file.path(extdataDir, "bcbio")

bcb <- loadRNASeq(
    uploadDir = uploadDir,
    interestingGroups = "group",
    ensemblVersion = 90)
metrics <- metrics(bcb)

dds <- bcbio(bcb, "DESeqDataSet")
design(dds) <- formula(~group)

rld <- rlog(dds)
res <- results(dds)

saveData(
    bcb,
    metrics,
    dds,
    rld,
    res,
    dir = extdataDir,
    compress = "xz",
    overwrite = TRUE)
