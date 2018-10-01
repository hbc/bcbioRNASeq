.meltCounts <- function(object, normalized) {
    assert_is_all_of(object, "bcbioRNASeq")
    assert_is_a_string(normalized)

    # Generate a SummarizedExperiment subset with only non-zero counts.
    # We're coercing to RSE here to enable fast subsetting on the object.
    rse <- as(object, "RangedSummarizedExperiment")
    # Note that we're comparing against the raw counts.
    keep <- rowSums(counts(rse)) > 0L
    rse <- rse[keep, , drop = FALSE]
    message(paste(
        nrow(rse), "/", nrow(object),
        "non-zero genes detected."
    ))

    # Subset the desired normalized counts, and slot into our SE object.
    counts <- counts(object, normalized = normalized) %>%
        .[keep, , drop = FALSE]
    assays(rse) <- list(counts = counts)

    # Apply scale transformation to counts axis, if necessary.
    trans <- .normalizedTrans(normalized)

    # Using our SummarizedExperiment method to return melted counts.
    data <- meltCounts(object = rse, trans = trans)
    assert_is_all_of(data, "grouped_df")

    data
}



# Determine if we should apply a log transformation to our normalized counts
# when plotting. Note that DESeqTransform are already log2 scale.
.normalizedTrans <- function(normalized) {
    if (normalized %in% c("rlog", "vst")) {
        # Already log2 scale.
        trans <- "identity"
    } else {
        trans <- "log2"
    }
}
