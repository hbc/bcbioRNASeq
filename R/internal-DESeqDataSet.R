.regenerateDESeqDataSet <- function(rse) {
    assert_is_all_of(rse, "RangedSummarizedExperiment")
    message(paste(
        "Generating DESeqDataSet with DESeq2",
        packageVersion("DESeq2")
    ))
    txi <- .regenerateTximportList(rse)
    dds <- DESeqDataSetFromTximport(
        txi = txi,
        colData = colData(rse),
        # Use an empty design formula
        design = ~ 1  # nolint
    )
    # Suppress warning about empty design formula
    dds <- suppressWarnings(DESeq(dds))
    validObject(dds)
    dds
}
