#' `DESeqResultsTables` Generator
#'
#' @section Obtaining results from DESeq2:
#'
#' It is recommended to specify the `contrast` argument as `character`:
#'
#' 1. Design matrix factor of interest.
#' 2. Numerator for LFC (expt).
#' 3. Denominator for LFC (control).
#'
#' For example, to get the relative change in gene expression between mutant
#' and wildtype samples, here's how to set the contrast vector:
#'
#' ```
#' c(
#'     factor = "genotype",
#'     numerator = "mutant",
#'     denominator = "wildtype"
#' )
#' ```
#'
#' @section Log fold change threshold cutoffs:
#'
#' Refer to [DESeq2::results()] for additional information about using
#' `lfcThreshold` and `altHypothesis` to set an alternative hypothesis based on
#' expected fold changes. In addition, the "Hypothesis tests with thresholds on
#' effect size" section in the DESeq2 paper provides additional explanation.
#'
#' Don't apply *post hoc* LFC threshold filtering to obtain results tables, as
#' this will destroy the meaning of the adjusted *P* values computed. If you are
#' expecting a biological effect size past a particular threshold, set
#' `lfcThreshold`, but be conservative.
#'
#' This [thread][] on the Bioconductor forums explains how [DESeq2::results()]
#' should be called with regard to LFC cutoffs in nice detail. In particular,
#' refer to Mike Love:
#'
#' "A common procedure is to disregard genes whose estimated LFC *β ir* is below
#' some threshold, *β ir ≤ θ*. However, this approach loses the benefit of an
#' easily interpretable FDR, as the reported *P* value and adjusted *P* value
#' still correspond to the test of *zero* LFC. It is therefore desirable to
#' include the threshold in the statistical testing procedure directly, i.e.,
#' not to filter post hoc on a reported fold-change *estimate*, but rather to
#' evaluate statistically directly whether there is sufficient evidence that the
#' LFC is above the chosen threshold."
#'
#' [thread]: https://support.bioconductor.org/p/101504/
#'
#' @name DESeqResultsTables
#' @family S4 Generators
#' @author Michael Steinbaugh
#' @include AllGenerics.R
#' @export
#'
#' @inheritParams general
#' @param rowData `boolean`. Include gene annotations.
#' @param counts `boolean`. Include DESeq2 normalized counts.
#'
#' @return `DESeqResultsTables`.
#'
#' @seealso
#' - [DESeq2::results()].
#' - [markdown()], [write()].
#'
#' @examples
#' data(deseq_small)
#'
#' # DESeqAnalysis ====
#' # This is the recommended default method.
#' x <- DESeqResultsTables(deseq_small)
#' print(x)
#'
#' # DESeqResults ====
#' res <- as(deseq_small, "DESeqResults")
#' x <- DESeqResultsTables(res)
#' print(x)
NULL



.DESeqResultsTables.DESeqResults <-  # nolint
    function(object) {
        validObject(object)
        assert_is_all_of(object, "DESeqResults")
        assert_is_subset(c("log2FoldChange", "padj"), colnames(object))
        alpha <- metadata(object)[["alpha"]]
        assertIsAlpha(alpha)
        lfcThreshold <- metadata(object)[["lfcThreshold"]]
        assert_is_a_number(lfcThreshold)
        assert_all_are_non_negative(lfcThreshold)

        # Set LFC and test (P value) columns.
        lfcCol <- "log2FoldChange"
        testCol <- "padj"
        lfc <- sym(lfcCol)
        test <- sym(testCol)
        assert_is_subset(
            x = c(lfcCol, testCol),
            y = colnames(object)
        )

        # DEG tables are sorted by adjusted P value.
        deg <- object %>%
            as_tibble(rownames = "rowname") %>%
            # Remove genes without an adjusted P value.
            filter(!is.na(!!test)) %>%
            # Remove genes that don't pass alpha cutoff.
            filter(!!test < !!alpha) %>%
            # Arrange by adjusted P value.
            arrange(!!test) %>%
            # Remove genes that don't pass LFC threshold.
            filter(!!lfc > !!lfcThreshold | !!lfc < -UQ(lfcThreshold))
        # Get directional subsets.
        degUp <- filter(deg, !!lfc > 0L)
        degDown <- filter(deg, !!lfc < 0L)

        new(
            Class = "DESeqResultsTables",
            all = object,
            deg = as(deg, "DataFrame"),
            degUp = as(degUp, "DataFrame"),
            degDown = as(degDown, "DataFrame")
        )
    }



.DESeqResultsTables.DESeqAnalysis <-  # nolint
    function(
        object,
        results = 1L,
        lfcShrink = TRUE,
        rowData = TRUE,
        counts = TRUE
    ) {
        validObject(object)
        assert_is_a_bool(rowData)
        assert_is_a_bool(counts)

        # Match the DESeqResults object.
        results <- .matchResults(
            object = object,
            results = results,
            lfcShrink = lfcShrink
        )

        # Add columns (optional) -----------------------------------------------
        if (isTRUE(rowData) || isTRUE(counts)) {
            # Coerce to `DataFrame`.
            # We'll regenerate a modified `DESeqResults` from this below.
            data <- as(results, "DataFrame")

            # Row annotations.
            if (isTRUE(rowData)) {
                message(paste(
                    "Adding `rowData()` annotations (atomic columns only)."
                ))
                rowData <- sanitizeRowData(rowData(object@data))
                # DESeq2 includes additional information in `rowData()` that
                # isn't informative for a user, and doesn't need to be included
                # in the CSV. Use our `bcb_small` example dataset to figure out
                # which columns are worth including.
                keep <- intersect(
                    x = colnames(rowData),
                    y = colnames(rowData(bcbioRNASeq::bcb_small))
                )
                rowData <- rowData[, keep, drop = FALSE]
                assert_all_are_true(vapply(rowData, is.atomic, logical(1L)))
                assert_is_non_empty(rowData)
                assert_are_disjoint_sets(colnames(data), colnames(rowData))
                data <- cbind(data, rowData)
            }

            # DESeq2 normalized counts
            if (isTRUE(counts)) {
                message("Adding DESeq2 normalized counts.")
                counts <- counts(object@data, normalized = TRUE)
                assert_are_disjoint_sets(colnames(data), colnames(counts))
                assert_are_identical(rownames(data), rownames(counts))
                data <- cbind(data, counts)
            }

            # Regenerate the DESeqResults.
            results <- DESeqResults(
                DataFrame = data,
                priorInfo = priorInfo(results)
            )
        }

        # Return ---------------------------------------------------------------
        DESeqResultsTables(results)
    }



#' @rdname DESeqResultsTables
#' @export
setMethod(
    f = "DESeqResultsTables",
    signature = signature("DESeqResults"),
    definition = .DESeqResultsTables.DESeqResults
)



#' @rdname DESeqResultsTables
#' @export
setMethod(
    f = "DESeqResultsTables",
    signature = signature("DESeqAnalysis"),
    definition = .DESeqResultsTables.DESeqAnalysis
)
