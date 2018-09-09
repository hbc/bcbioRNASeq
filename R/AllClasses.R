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
        all = "DESeqResults",
        deg = "DataFrame",
        degUp = "DataFrame",
        degDown = "DataFrame"
    )
)

setValidity(
    Class = "DESeqResultsTables",
    method = function(object) {
        assert_is_a_string(contrastName(object@all))
        assert_is_a_number(metadata(object@all)[["alpha"]])
        assert_is_a_number(metadata(object@all)[["lfcThreshold"]])
        assert_is_subset(
            x = rownames(object@downregulated),
            y = rownames(object@all)
        )
        assert_is_subset(
            x = rownames(object@upregulated),
            y = rownames(object@all)
        )
        assert_are_disjoint_sets(
            x = rownames(object@upregulated),
            y = rownames(object@downregulated)
        )

        # Require that the subsets match the `DESeqResults` summary.
        match <- removeNA(str_match(
            string = capture.output(summary(object@all)),
            pattern = "^LFC.*\\s\\:\\s([0-9]+).*"
        ))
        assert_are_identical(
            x = nrow(object@upregulated),
            y = as.integer(match[1, 2])
        )
        assert_are_identical(
            x = nrow(object@downregulated),
            y = as.integer(match[2, 2])
        )

        TRUE
    }
)
