#' Count Matrix Accessors
#'
#' By default, [counts()] returns the raw counts. Normalized counts, including
#' transcripts per million (TPM) can be accessed using the `normalized`
#' argument.
#'
#' @name counts
#' @family Data Functions
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
#' @return `matrix`.
#'
#' @examples
#' # Raw counts (from tximport)
#' counts(bcb_small, normalized = FALSE) %>% summary()
#'
#' # TPM (from tximport)
#' counts(bcb_small, normalized = "tpm") %>% summary()
#'
#' # Normalized counts (from DESeq2)
#' counts(bcb_small, normalized = TRUE) %>% summary()
#'
#' # rlog (from DESeq2)
#' counts(bcb_small, normalized = "rlog") %>% summary()
#'
#' # VST (from DESeq2)
#' counts(bcb_small, normalized = "vst") %>% summary()
#'
#' # TMM (from edgeR)
#' counts(bcb_small, normalized = "tmm") %>% summary()
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
            if (identical(normalized, FALSE)) {
                # Raw counts from tximport (default)
                counts <- assay(object)
            } else if (identical(normalized, TRUE)) {
                # DESeq2 normalized counts
                dds <- as(object, "DESeqDataSet")
                assert_is_all_of(dds, "DESeqDataSet")
                validObject(dds)
                counts <- counts(dds, normalized = TRUE)
            }
        } else if (is.character(normalized)) {
            assert_is_a_string(normalized)
            assert_is_subset(
                x = normalized,
                y = c("raw", "tpm", "rlog", "vst", "tmm")
            )
            if (normalized == "tmm") {
                # Calculate TMM on the fly
                counts <- tmm(assay(object))
            } else {
                # Use matrices slotted into `assays()`
                counts <- assays(object)[[normalized]]
            }
            if (!is.matrix(counts)) {
                warn(paste(
                    normalized, "not present in assays.",
                    "Calculating log2 TMM counts instead."
                ))
                counts <- tmm(assay(object))
                inform("Applying log2 transformation to TMM values")
                counts <- log2(counts + 1L)
            }
        }

        assert_is_matrix(counts)
        counts
    }
)
