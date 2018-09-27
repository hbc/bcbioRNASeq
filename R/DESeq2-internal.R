# Calculate a numeric vector to define the colors
# -1: downregulated
#  0: not significant
#  1: upregulated
.addIsDECol <- function(
    data,
    testCol = "padj",
    alpha,
    lfcCol = "log2FoldChange",
    lfcThreshold = 0L
) {
    # test: P value or S value
    test <- data[[testCol]]
    # lfc: log2 fold change cutoff
    lfc <- data[[lfcCol]]
    isDE <- mapply(
        test = test,
        lfc = lfc,
        FUN = function(test, lfc) {
            if (any(is.na(c(test, lfc)))) {
                # nonsignificant
                0L
            } else if (test < alpha & lfc > lfcThreshold) {
                # upregulated
                1L
            } else if (test < alpha & lfc < -lfcThreshold) {
                # downregulated
                -1L
            } else {
                0L
            }
        },
        SIMPLIFY = TRUE,
        USE.NAMES = FALSE
    )
    isDE <- as.factor(isDE)
    data[["isDE"]] <- isDE
    data
}



.matchResults <- function(
    object,
    results,
    lfcShrink = FALSE
) {
    assert_is_all_of(object, "DESeqAnalysis")
    # Default to using the first contrast, for convenience.
    if (missing(results)) {
        results <- 1L
    }
    assert_is_scalar(results)
    assert_is_a_bool(lfcShrink)
    if (isTRUE(lfcShrink)) {
        slotName <- "lfcShrink"
    } else {
        slotName <- "results"
    }
    results <- slot(object, name = slotName)[[results]]
    assert_is_all_of(results, "DESeqResults")

    # Inform the user about which data we're using.
    msg <- paste(bold("DESeqResults:"), contrastName(results))
    if (isTRUE(lfcShrink)) {
        msg <- paste(msg, "(shrunken LFC)")
    }
    message(msg)

    results
}



.regenerateDESeqDataSet <- function(object) {
    assert_is_all_of(object, "RangedSummarizedExperiment")
    message(paste0(
        "Generating DESeqDataSet with DESeq2 ",
        packageVersion("DESeq2"), "..."
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
