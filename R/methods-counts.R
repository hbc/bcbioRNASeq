#' Count Matrix Accessors
#'
#' By default, [counts()] returns the raw counts. Normalized counts, including
#' transcripts per million (TPM) can be accessed using the "`normalized`"
#' argument.
#'
#' @name counts
#' @family Data Functions
#' @author Michael Steinbaugh, Lorena Pantano
#'
#' @importFrom BiocGenerics counts
#'
#' @inheritParams general
#' @param normalized Logical or character indicating which normalization
#'   method to apply:
#'   - `FALSE`: Raw counts (tximport).
#'   - `TRUE`: DESeq2 normalized counts. Calculated on the fly.
#'   - "`tpm`": Transcripts per million (tximport).
#'   - "`tmm`": edgeR trimmed mean of M-values. Calculated on the fly.
#'   - "`rlog`": DESeq2 **log2** regularized log transformation.
#'   - "`vst`": DESeq2 **log2** variance stabilizing transformation.
#'
#' @return `matrix`.
#'
#' @seealso
#' - [tpm()].
#' - [tmm()].
#' - [DESeq2::counts()].
#' - [DESeq2::rlog()].
#' - [DESeq2::varianceStabilizingTransformation()].
#'
#' @examples
#' counts(bcb_small, normalized = FALSE) %>% summary()
#' counts(bcb_small, normalized = TRUE) %>% summary()
#' counts(bcb_small, normalized = "tpm") %>% summary()
#' counts(bcb_small, normalized = "tmm") %>% summary()
#'
#' # log2
#' counts(bcb_small, normalized = "rlog") %>% summary()
#' counts(bcb_small, normalized = "vst") %>% summary()
NULL



# Methods ======================================================================
#' @rdname counts
#' @export
setMethod(
    "counts",
    signature("bcbioRNASeq"),
    function(object, normalized = FALSE) {
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
                # Support for skipped DESeq2 transforms: log2 TMM
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
