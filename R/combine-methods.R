# Here we are using the combine method for SummarizedExperiment defined in the
# basejump package. This method is strict, and requires that the `metadata()`
# slot for the two objects is identical. So here we must first prepare the
# bcbioRNASeq objects with more minimal metadata prior to the `combine()` call.



#' @name combine
#' @inherit basejump::combine
#' @importMethodsFrom basejump combine
#' @return `RangedSummarizedExperiment`.
NULL



#' @rdname combine
#' @name combine
#' @importFrom BiocGenerics combine
#' @usage combine(x, y, ...)
#' @export
NULL



combine.bcbioRNASeq <-  # nolint
    function(x, y) {
        validObject(x)
        validObject(y)

        message(paste(
            "Combining bcbioRNASeq runs.",
            "Note that this step returns as RangedSummarizedExperiment.",
            sep = "\n"
        ))

        # Coerce objects to RSE so we can modify `metadata()` slot.
        x <- as(x, "RangedSummarizedExperiment")
        y <- as(y, "RangedSummarizedExperiment")

        # Metadata -------------------------------------------------------------
        message("Updating metadata.")

        # Check that the metadata names (but not the content) are identical.
        assert(identical(
            x = names(metadata(x)),
            y = names(metadata(y))
        ))

        mx <- metadata(x)
        my <- metadata(y)

        # Strip out attributes on dataVersions prior to comparison, otherwise
        # these won't return identical, due to differing brio import attributes.
        assert(
            is(mx[["dataVersions"]], "DataFrame"),
            is(my[["dataVersions"]], "DataFrame")
        )
        attr(mx[["dataVersions"]]@listData, "brio") <- NULL
        attr(my[["dataVersions"]]@listData, "brio") <- NULL

        # Drop blacklisted metadata elements that are expected to differ between
        # the objects, and therefore shouldn't be kept in the merge.
        blacklist <- c(
            "allSamples",
            "bcbioCommandsLog",
            "bcbioLog",
            "call",
            "date",
            "projectDir",
            "runDate",
            "sampleDirs",
            "sampleMetadataFile",
            "sessionInfo",
            "uploadDir",
            "wd",
            "yaml"
        )

        mx <- mx[setdiff(names(mx), blacklist)]
        my <- my[setdiff(names(my), blacklist)]
        metadata(x) <- mx
        metadata(y) <- my
        assert(identical(metadata(x), metadata(y)))

        # Column data ----------------------------------------------------------
        message("Updating column data.")

        # Check for column mismatches and restore NA values, if necessary. This
        # mismatch can occur because our metadata importer will drop columns
        # with all NA values, which is useful for handling human metadata. This
        # can create a column mismatch when we're subsetting large sequencing
        # runs into batches.

        cdx <- colData(x)
        cdy <- colData(y)

        cols <- union(names(cdx), names(cdy))

        if (!isTRUE(setequal(cols, names(cdx)))) {
            diff <- setdiff(cols, names(cdx))
            for (col in diff) {
                cdx[[col]] <- NA
            }
        }

        if (!isTRUE(setequal(cols, names(cdy)))) {
            diff <- setdiff(cols, names(cdy))
            for (col in diff) {
                cdy[[col]] <- NA
            }
        }

        colData(x) <- cdx
        colData(y) <- cdy

        # Return ---------------------------------------------------------------
        combine(x, y)
    }



#' @rdname combine
#' @export
setMethod(
    f = "combine",
    signature = signature(
        x = "bcbioRNASeq",
        y = "bcbioRNASeq"
    ),
    definition = combine.bcbioRNASeq
)
