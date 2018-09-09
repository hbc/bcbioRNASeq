setValidity(
    "DESeqAnalysis",
    function(object) {
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
