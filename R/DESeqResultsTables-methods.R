#' @name DESeqResultsTables
#' @inherit DESeqResultsTables-class
#' @inheritParams params
#' @inheritParams basejump::params
#' @author Michael Steinbaugh
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
#' "A common procedure is to disregard genes whose estimated LFC *beta ir* is
#' below some threshold, *beta ir ≤ theta*. However, this approach loses the
#' benefit of an easily interpretable FDR, as the reported *P* value and
#' adjusted *P* value still correspond to the test of *zero* LFC. It is
#' therefore desirable to include the threshold in the statistical testing
#' procedure directly, i.e., not to filter post hoc on a reported fold-change
#' *estimate*, but rather to evaluate statistically directly whether there is
#' sufficient evidence that the LFC is above the chosen threshold."
#'
#' [thread]: https://support.bioconductor.org/p/101504/
#'
#' @return `DESeqResultsTables`.
#'
#' @seealso
#' - [DESeq2::results()].
#' - [markdown()], [write()].
#'
#' @examples
#' data(deseq)
#'
#' ## DESeqAnalysis ====
#' ## This is the recommended default method.
#' x <- DESeqResultsTables(deseq)
#' print(x)
#'
#' ## DESeqResults ====
#' res <- as(deseq, "DESeqResults")
#' x <- DESeqResultsTables(res)
#' print(x)
NULL



# DESeqResults =================================================================
DESeqResultsTables.DESeqResults <-  # nolint
    function(object) {
        validObject(object)
        assert_that(is(object, "DESeqResults"))
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
        up <- deg %>%
            filter(!!lfc > 0L) %>%
            pull("rowname")
        down <- deg %>%
            filter(!!lfc < 0L) %>%
            pull("rowname")

        new(
            Class = "DESeqResultsTables",
            results = object,
            deg = list(up = up, down = down)
        )
    }

# DESeqAnalysis ================================================================
DESeqResultsTables.DESeqAnalysis <-  # nolint
    function(
        object,
        results = 1L,
        lfcShrink = TRUE
    ) {
        validObject(object)

        # Prepare the DESeqResultsTables object with our DESeqResults method.
        results <- .matchResults(
            object = object,
            results = results,
            lfcShrink = lfcShrink
        )
        # Prepare the DESeqResultsTables object.
        out <- DESeqResultsTables(results)
        rm(results)

        # Automatically populate additional slots using DESeqDataSet.
        data <- as(object, "DESeqDataSet")

        # Note that we're slotting the size factor-normalized counts here.
        counts <- counts(data, normalized = TRUE)
        slot(out, "counts") <- counts

        # Slot the row annotations (genomic ranges).
        #
        # DESeq2 includes additional columns in `mcols()` that aren't
        # informative for a user, and doesn't need to be included in the tables.
        #
        # Instead, only keep informative columns that are character or factor.
        # Be sure to drop complex, non-atomic columns (e.g. list, S4) that are
        # allowed in GRanges/DataFrame but will fail to write to disk as CSV.
        rowRanges <- rowRanges(data)
        mcols <- mcols(rowRanges)
        keep <- vapply(
            X = mcols,
            FUN = function(x) {
                is.character(x) || is.factor(x)
            },
            FUN.VALUE = logical(1L)
        )
        mcols <- mcols[, keep, drop = FALSE]
        mcols(rowRanges) <- mcols
        assert_is_non_empty(rowRanges)
        assert_are_identical(
            x = rownames(data),
            y = names(rowRanges)
        )
        assert_are_disjoint_sets(
            x = colnames(data),
            y = colnames(mcols(rowRanges))
        )
        assert_all_are_true(vapply(
            X = mcols(rowRanges),
            FUN = is.atomic,
            FUN.VALUE = logical(1L)
        ))
        slot(out, "rowRanges") <- rowRanges

        # Slot human-friendly sample names, if they are defined.
        sampleNames <- sampleNames(data)
        if (
            has_length(sampleNames) &&
            !identical(
                x = as.character(sampleNames),
                y = colnames(data)
            )
        ) {
            slot(out, "sampleNames") <- sampleNames
        }


        # Slot metadata.
        slot(out, "metadata") <- list(
            version = metadata(data)[["version"]]
        )

        out
    }



#' @rdname DESeqResultsTables
#' @export
setMethod(
    f = "DESeqResultsTables",
    signature = signature("DESeqResults"),
    definition = DESeqResultsTables.DESeqResults
)



#' @rdname DESeqResultsTables
#' @export
setMethod(
    f = "DESeqResultsTables",
    signature = signature("DESeqAnalysis"),
    definition = DESeqResultsTables.DESeqAnalysis
)
