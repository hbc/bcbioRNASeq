devtools::load_all()

extdataDir <- file.path("inst", "extdata")

bcb <- loadRNASeq(
    uploadDir = file.path(extdataDir, "bcbio"),
    interestingGroups = "group",
    ensemblVersion = 90)

dds <- bcbio(bcb, "DESeqDataSet")
design(dds) <- formula(~group)
dds <- DESeq(dds)

rld <- rlog(dds)
res <- results(dds, contrast = c("group", "ko", "ctrl"))

saveData(
    bcb,
    dds,
    rld,
    res,
    dir = extdataDir,
    compress = "xz",
    overwrite = TRUE)
