setClassUnion(
    name = "missingOrNULL",
    members = c("missing", "NULL")
)



# bcbioRNASeq ==================================================================
#' @rdname bcbioRNASeq
#' @aliases NULL
#' @exportClass bcbioRNASeq
#' @usage NULL
setClass(
    Class = "bcbioRNASeq",
    contains = "RangedSummarizedExperiment"
)



# DESeqAnalysis ================================================================
#' @rdname DESeqAnalysis
#' @aliases NULL
#' @exportClass DESeqAnalysis
#' @usage NULL
setClass(
    Class = "DESeqAnalysis",
    slots = list(
        DESeqDataSet = "DESeqDataSet",
        DESeqTransform = "DESeqTransform",
        DESeqResults = "list",
        lfcShrink = "list"
    )
)

setValidity(
    Class = "DESeqAnalysis",
    method = function(object) {
        # Require valid dimnames for data set.
        assertHasValidDimnames(object@DESeqDataSet)

        # Require gene-to-symbol mappings.
        assert_is_all_of(
            x = gene2symbol(object@DESeqDataSet),
            classes = "gene2symbol"
        )

        # DESeqDataSet and DESeqTransform must match.
        assert_are_identical(
            x = dimnames(object@DESeqDataSet),
            y = dimnames(object@DESeqTransform)
        )

        # DESeqDataSet and DESeqResults must match.
        invisible(lapply(
            X = object@DESeqResults,
            FUN = function(res) {
                assert_are_identical(
                    x = rownames(res),
                    y = rownames(object@DESeqDataSet)
                )
            }
        ))

        # DESeqResults and lfcShrink rownames must match.
        invisible(mapply(
            unshrunken = object@DESeqResults,
            shrunken = object@lfcShrink,
            FUN = function(unshrunken, shrunken) {
                assert_are_identical(
                    x = rownames(unshrunken),
                    y = rownames(shrunken)
                )
            }
        ))

        TRUE
    }
)



# resultsTables ================================================================
#' @rdname resultsTables
#' @aliases NULL
#' @exportClass resultsTables
#' @usage NULL
setClass(
    Class = "DESeqResultsTables",
    slots = list(
        deg = "tbl_df",
        degLFC = "tbl_df",
        degLFCUp = "tbl_df",
        degLFCDown = "tbl_df",
        all = "DESeqResults",
        contrast = "character",
        alpha = "numeric",
        lfcThreshold = "numeric"
    )
)

# FIXME Improve this.
setValidity(
    Class = "DESeqResultsTables",
    method = function(object) {
        assert_is_subset(
            x = c("all", "deg", "degLFCDown", "degLFCUp"),
            y = names(object)
        )
        TRUE
    }
)
