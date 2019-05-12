#' Count matrix accessors
#'
#' By default, [counts()] returns the raw counts. Normalized counts, including
#' transcripts per million (TPM) can be accessed using the "`normalized`"
#' argument.
#'
#' @name counts
#' @family Data Functions
#' @author Michael Steinbaugh, Lorena Pantano
#'
#' @inheritParams general
#' @param normalized `string` or `boolean`. Which normalization method to apply:
#'   - `FALSE`: Raw counts (tximport).
#'   - `TRUE`: DESeq2 normalized counts. Calculated on the fly.
#'   - "`tpm`": Transcripts per million (tximport).
#'   - "`vst`": DESeq2 **log2** variance stabilizing transformation.
#'   - "`rlog`": DESeq2 **log2** regularized log transformation.
#'   - "`tmm`": edgeR trimmed mean of M-values. Calculated on the fly.
#'   - "`rle`": Relative log expression transformation.
#' @param ... Additional arguments.
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
#' # bcbioRNASeq ====
#' counts(bcb_small, normalized = FALSE) %>% summary()
#' counts(bcb_small, normalized = TRUE) %>% summary()
#' counts(bcb_small, normalized = "tpm") %>% summary()
#' counts(bcb_small, normalized = "vst") %>% summary()
NULL



#' @rdname counts
#' @name counts
#' @importFrom BiocGenerics counts
#' @usage counts(object, ...)
#' @export
NULL



counts.bcbioRNASeq <-  # nolint
    function(object, normalized = FALSE) {
        validObject(object)
        assert_is_any_of(normalized, c("character", "logical"))

        if (is.logical(normalized)) {
            if (identical(normalized, FALSE)) {
                counts <- assay(object)
            } else if (identical(normalized, TRUE)) {
                counts <- assays(object)[["normalized"]]
            }
        } else if (is.character(normalized)) {
            assert_is_a_string(normalized)
            assert_is_subset(
                x = normalized,
                y = c("tpm", "vst", "rlog", "tmm", "rle")
            )
            if (normalized == "tmm") {
                # Calculate TMM on the fly
                counts <- tmm(assay(object))
            } else if (normalized == "rle") {
                # Calculate RLE on the fly
                counts <- t(t(assay(object)) / colMedians(assay(object)))
            } else {
                # Use matrices slotted into `assays()`
                counts <- assays(object)[[normalized]]
            }
            if (!is.matrix(counts)) {
                # Support for skipped DESeq2 transforms: log2 TMM
                warning(paste(
                    normalized, "not present in assays.",
                    "Calculating log2 TMM counts instead."
                ), call. = FALSE)
                counts <- tmm(assay(object))
                message("Applying log2 transformation to TMM values")
                counts <- log2(counts + 1L)
            }
        }

        assert_is_matrix(counts)
        counts
    }



#' @rdname counts
#' @export
setMethod(
    f = "counts",
    signature = signature("bcbioRNASeq"),
    definition = counts.bcbioRNASeq
)
