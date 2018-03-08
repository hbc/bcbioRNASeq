#' Count Matrix Accessors
#'
#' By default, [counts()] returns the raw counts. Normalized counts, including
#' transcripts per million (TPM) can be accessed using the `normalized`
#' argument.
#'
#' @rdname counts
#' @name counts
#' @author Michael Steinbaugh
#'
#' @importFrom BiocGenerics counts
#'
#' @inheritParams general
#'
#' @param normalized Select raw counts (`FALSE`), DESeq2 normalized counts
#'   (`TRUE`), or additional normalization methods:
#'   - `tpm`: Transcripts per million.
#'   - `tmm`: Trimmed mean of M-values (edgeR).
#'   - `rlog`: Regularized log transformation ([DESeq2::rlog()]).
#'   - `vst`: Variance stabilizing transformation
#'     ([DESeq2::varianceStabilizingTransformation()]).
#'
#' @return [matrix].
#'
#' @examples
#' load(system.file("extdata/bcb.rda", package = "bcbioRNASeq"))
#'
#' # Raw counts (from tximport)
#' counts(bcb, normalized = FALSE) %>% summary()
#'
#' # TPM (from tximport)
#' counts(bcb, normalized = "tpm") %>% summary()
#'
#' # Normalized counts (from DESeq2)
#' counts(bcb, normalized = TRUE) %>% summary()
#'
#' # rlog (from DESeq2)
#' counts(bcb, normalized = "rlog") %>% summary()
#'
#' # TMM (from edgeR)
#' counts(bcb, normalized = "tmm") %>% summary()
#'
#' # VST (from DESeq2)
#' counts(bcb, normalized = "vst") %>% summary()
NULL



# Methods ======================================================================
#' @rdname counts
#' @export
setMethod(
    "counts",
    signature("bcbioRNASeq"),
    function(
        object,
        normalized = FALSE
) {
        validObject(object)
        assert_is_any_of(normalized, c("character", "logical"))

        if (is.logical(normalized)) {
            if (isTRUE(normalized)) {
                # DESeq2 normalized counts
                dds <- assays(object)[["dds"]]
                assert_is_all_of(dds, "DESeqDataSet")
                validObject(dds)
                counts <- counts(dds, normalized = TRUE)
            } else {
                # Raw counts from tximport
                counts <- assays(object)[["raw"]]
            }
        } else if (is.character(normalized)) {
            assert_is_a_string(normalized)
            assert_is_subset(
                normalized,
                c(
                    "tpm",
                    "tmm",
                    "rlog",
                    "vst"
                )
            )
            counts <- assays(object)[[normalized]]

            # Return the matrix from slotted DESeqTransform
            if (is(counts, "DESeqTransform")) {
                counts <- assay(counts)
            }

            # Return log2 TMM if assay is NULL.
            # This helps support skipping CPU-intensive DESeq2 variance
            # stabilization calculations on large datasets.
            if (is.null(counts)) {
                warn(paste(
                    normalized, "counts not defined.",
                    "Calculating and using log2 tmm counts on the fly instead."
                ))
                counts <- tmm(object)
                counts <- log2(counts + 1L)
            }
        }

        assert_is_matrix(counts)
        counts
    }
)
