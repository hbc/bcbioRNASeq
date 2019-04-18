.dynamicTrans <- function(
    object,
    normalized,
    trans
) {
    validObject(object)
    normalized <- match.arg(normalized, choices = normalizedCounts)
    trans <- match.arg(trans, choices = c("log2", "log10"))

    # Don't allow DESeqTransform to be plotted on log10 scale.
    if (
        normalized %in% c("vst", "rlog") &&
        trans != "log2"
    ) {
        message(paste(normalized, "counts are log2 scale."))
    }

    # Check for DESeqTransform that are already log2 scale and update `trans`.
    if (normalized %in% c("rlog", "vst")) {
        trans <- "identity"
    }

    # Coerce to RangedSummarizedExperiment.
    counts <- counts(object, normalized = normalized)
    rse <- as(object, "RangedSummarizedExperiment")
    assays(rse) <- list(counts)
    assayNames(rse) <- normalized

    # Set the counts axis label.
    countsAxisLabel <- paste(normalized, "counts")
    if (trans != "identity") {
        countsAxisLabel <- paste(trans, countsAxisLabel)
    }

    list(
        object = rse,
        trans = trans,
        countsAxisLabel = countsAxisLabel
    )
}



# nolint start
#
# Reverse the order of a categorical axis in ggplot2.
# - https://gist.github.com/jennybc/6f3fa527b915b920fdd5
# - https://gist.github.com/jennybc/6f3fa527b915b920fdd5#gistcomment-2709675
#
# Useful alternate approach:
# > aes(x = reorder(the_factor, desc(the_factor)), ...)
#
# nolint end

# Make the samples human readable when flipped onto Y axis.
.flipMode <- function(object) {
    assert(is(object, "ggplot"))
    data <- object[["data"]]
    assert(is.data.frame(data))
    samples <- data[["sampleName"]]
    assert(is.factor(samples))
    object +
        scale_x_discrete(limits = rev(levels(samples))) +
        coord_flip()
}
