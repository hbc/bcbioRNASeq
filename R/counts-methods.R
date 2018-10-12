# FIXME Consider dropping `tmm()` generic.



#' Count Matrix Accessors
#'
#' By default, returns the raw counts. Normalized counts in a variety of formats
#' can be accessed using the `normalized` argument.
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
#'   - `FALSE`: Raw counts.
#'     - When using a [tximport][]-compatible caller, these are length scaled
#'       by default (see `countsFromAbundance`).
#'     - When using a [featureCounts][]-compatible caller, these are `integer`.
#'
#' [tximport][] caller-specific normalizations:
#'
#'   - `"tpm"`: **T**ranscripts **p**er **m**illion.
#'
#' Additional gene-level-specific normalizations:
#'
#'   - `TRUE`: Size factor-adjusted counts.\cr
#'     See [DESeq2::sizeFactors] for more information.
#'   - `"vst"`: **V**ariance-**s**tabilizing **t**ransformation.\cr
#'     Requires `vst = TRUE` to be set during [bcbioRNASeq()] call.\cr
#'     See [DESeq2::varianceStabilizingTransformation] for more information.
#'   - `"rlog"`: **R**egularized **log** transformation (log2).\cr
#'     Requires `rlog = TRUE` to be set during [bcbioRNASeq()] call.\cr
#'     See [DESeq2::rlog] for more information.
#'   - `"tmm"`: **T**rimmed **m**ean of **M**-values.\cr
#'     See [edgeR::calcNormFactors] for more information.
#'   - `"rle"`: **R**elative **l**og **e**xpression transformation.
#'   - `"fpkm"`: **F**ragments **p**er **k**ilobase per **m**illion mapped
#'     fragments.\cr
#'     Requires annotations in [rowRanges()] with defined [widths], otherwise
#'     will be skipped during the [bcbioRNASeq()] load call.\cr
#'     See [DESeq2::fpkm] for more information.
#'
#' [featureCounts]: http://bioinf.wehi.edu.au/featureCounts/
#' [tximport]: https://bioconductor.org/packages/release/bioc/html/tximport.html
#'
#' @return `matrix`.
#'
#' @references
#' - TMM: Robinson and Oshlack (2010).
#' - RLE: Anders and Huber (2010).
#'
#' @seealso
#' - [tximport::tximport()].
#' - [DESeq2::counts()].
#' - [DESeq2::sizeFactors()].
#' - [DESeq2::varianceStabilizingTransformation()].
#' - [DESeq2::rlog()].
#' - [DESeq2::fpkm()].
#' - [edgeR::calcNormFactors()].
#' - [tmm()], [relativeLogExpression()].
#'
#' @examples
#' data(bcb_small)
#' counts(bcb_small) %>% summary()
NULL



.counts.bcbioRNASeq <-  # nolint
    function(object, normalized = FALSE) {
        validObject(object)
        assert_is_any_of(normalized, c("character", "logical"))
        # Ensure that primary assay matches counts.
        assert_are_identical(assayNames(object)[[1L]], "counts")

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
            counts <- assays(object)[["counts"]]
            counts <- tmm(counts)
        } else if (assayName == "rle") {
            counts <- assays(object)[["counts"]]
            counts <- relativeLogExpression(counts)
        } else {
            # Get matrix slotted in `assays()`.
            # Note that we've killed the log2 TMM fall back support if
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
