devtools::load_all()
extdata_dir <- file.path("inst", "extdata")

loadRemoteData("http://bcbiornaseq.seq.cloud/f1000v2/bcb_all.rda")

# Select only the control and day 7 folic acid samples
bcb <- selectSamples(bcb, day = c(0L, 7L))

# TODO Select only the top 500 most abundant genes

dds <- bcbio(bcb, "DESeqDataSet")
design(dds) <- formula(~day)
dds <- DESeq(dds)
rld <- rlog(dds)
res <- results(
    dds,
    contrast = c(
        factor = "day",
        numerator = 7L,
        denominator = 0L)
)

saveData(
    bcb,
    dds,
    rld,
    res,
    dir = extdata_dir,
    compress = "xz",
    overwrite = TRUE)
