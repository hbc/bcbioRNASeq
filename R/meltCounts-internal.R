.filterNonzeroGenes <- function(object) {
    assert_is_all_of(object, "bcbioRNASeq")
    keep <- rowSums(counts(object)) > 0L
    object <- object[keep, , drop = FALSE, transform = FALSE]
    message(paste(nrow(object), "non-zero genes detected"))
    object
}



.meltCounts <- function(object, normalized) {
    assert_is_all_of(object, "bcbioRNASeq")
    assert_is_a_string(normalized)

    # Apply log2 transformation, if  necessary.
    if (normalized %in% c("rlog", "vst")) {
        log2Transform <- FALSE
    } else {
        log2Transform <- TRUE
    }

    # Subset the object to only include non-zero genes.
    nonzero <- .filterNonzeroGenes(object)
    counts <- counts(nonzero, normalized = normalized)

    # Using SummarizedExperiment method to return melted counts.
    rse <- as(nonzero, "RangedSummarizedExperiment")
    assays(rse) <- list(counts = counts)
    data <- meltCounts(
        object = rse,
        log2Transform = log2Transform
    )
    assert_is_all_of(data, "grouped_df")
    data
}
