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
        upregulated = "DataFrame",
        downregulated = "DataFrame",
        all = "DESeqResults"
    )
)

setValidity(
    Class = "DESeqResultsTables",
    method = function(object) {
        # results <- slot(object, "DESeqResults")
        # assert_is_all_of(results, "DESeqResults")
        #
        # contrastName <- contrastName(results)
        # assert_is_a_string(contrastName)
        #
        # alpha <- metadata(results)[["alpha"]]
        # assert_is_a_number(alpha)
        #
        # lfcThreshold <- metadata(results)[["lfcThreshold"]]
        # assert_is_a_number(lfcThreshold)
        #
        # degUp <- slot(object, "degUp")
        # degDown <- slot(object, "degDown")
        # assert_is_subset(degUp, rownames(results))
        # assert_is_subset(degDown, rownames(results))
        # assert_are_disjoint_sets(degUp, degDown)

        TRUE
    }
)
