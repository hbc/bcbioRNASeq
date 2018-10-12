#' Count Matrix Accessors
#'
#' By default, [counts()] returns the raw counts. Normalized counts, including
#' transcripts per million (TPM) can be accessed using the `"normalized"`
#' argument.
#'
#' @name counts
#' @family Data Functions
#' @author Michael Steinbaugh, Lorena Pantano
#' @importFrom BiocGenerics counts counts<-
#' @export
#'
#' @inheritParams general
#' @param normalized `string` or `boolean`. Which normalization method to apply:
#'
#'   - `FALSE`: Raw counts (tximport).
#'   - `TRUE`: DESeq2 normalized counts. Calculated on the fly.
#'   - `"tpm"`: Transcripts per million (tximport).
#'   - `"vst"`: DESeq2 **log2** variance stabilizing transformation.
#'   - `"rlog"`: DESeq2 **log2** regularized log transformation.
#'   - `"tmm"`: edgeR trimmed mean of M-values. Calculated on the fly.
#'   - `"rle"`: Relative log expression transformation.
#'
#'   Note that transcript-level counts support only `FALSE` and `"tpm"`.
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
#' data(bcb_small)
#' counts(bcb_small, normalized = FALSE) %>% summary()
#' counts(bcb_small, normalized = TRUE) %>% summary()
#' counts(bcb_small, normalized = "tpm") %>% summary()
#' counts(bcb_small, normalized = "vst") %>% summary()
NULL



.counts.bcbioRNASeq <-  # nolint
    function(object, normalized = FALSE) {
        validObject(object)
        assert_is_any_of(normalized, c("character", "logical"))

        # Restrict the `normalized` arguments for transcript-level objects.
        if (.isTranscriptLevel(object)) {
            assert_is_subset(
                x = normalized,
                y = list(FALSE, "tpm")
            )
        }

        if (identical(normalized, FALSE)) {
            assayName <- "counts"
        } else if (identical(normalized, TRUE)) {
            assayName <- "normalized"
        } else {
            assayName <- match.arg(
                arg = normalized,
                choices = normalizedCounts
            )
        }
        assert_is_a_string(assayName)

        if (assayName == "tmm") {
            # Calculate TMM on the fly.
            counts <- tmm(assay(object))
        } else if (assayName == "rle") {
            # Calculate RLE on the fly.
            counts <- t(t(assay(object)) / colMedians(assay(object)))
        } else {
            # Get matrix slotted in `assays()`.
            # Note that we're killing the log2 TMM fall back support if
            # DESeq2 transforms are skipped, because that is confusing.
            counts <- assays(object)[[assayName]]
        }

        assert_is_matrix(counts)
        counts
    }



#' @rdname counts
#' @export
setMethod(
    f = "counts",
    signature = signature("bcbioRNASeq"),
    definition = .counts.bcbioRNASeq
)
