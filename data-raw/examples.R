devtools::load_all()

extdata_dir <- file.path("inst", "extdata")

bcb <- loadRNASeq(
    uploadDir = file.path(extdata_dir, "bcbio"),
    interestingGroups = "group",
    ensemblVersion = 90L)

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
    dir = extdata_dir,
    compress = "xz",
    overwrite = TRUE)
