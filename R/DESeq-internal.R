.regenerateDESeqDataSet <- function(object) {
    assert_is_all_of(object, "RangedSummarizedExperiment")
    message(paste(
        "Generating DESeqDataSet with DESeq2",
        packageVersion("DESeq2")
    ))
    txi <- .regenerateTximportList(object)
    dds <- DESeqDataSetFromTximport(
        txi = txi,
        colData = colData(object),
        # Use an empty design formula
        design = ~ 1L
    )
    # Suppress warning about empty design formula
    dds <- suppressWarnings(DESeq(dds))
    validObject(dds)
    dds
}



.transformCountsAxisLabel <- function(object) {
    paste(.transformType(object), "counts (log2)")
}



.transformType <- function(object) {
    assert_is_all_of(object, "DESeqTransform")
    if ("rlogIntercept" %in% colnames(mcols(object))) {
        "rlog"
    } else {
        # varianceStabilizingTransformation
        "vst"
    }
}
