#' Counts
#'
#' Count matrix.
#'
#' By default, return the raw counts. Normalized counts in a variety of formats
#' can be accessed using the `normalized` argument.
#'
#' @name counts
#' @author Michael Steinbaugh, Lorena Pantano
#' @note Updated 2019-08-07.
#'
#' @inheritParams acidroxygen::params
#' @param normalized `character(1)` or `logical(1)`.
#'   Which normalization method to apply:
#'
#'   - `FALSE`: Raw counts.
#'     When using a [tximport][]-compatible caller, these are length scaled
#'     by default (see `countsFromAbundance` argument).
#'     When using a [featureCounts][]-compatible caller, these are `integer`.
#'
#' [tximport][] caller-specific normalizations:
#'
#'   - `"tpm"`: **T**ranscripts **p**er **m**illion.
#'
#' Additional gene-level-specific normalizations:
#'
#'   - `TRUE` / `"sf"`: **S**ize **f**actor (i.e. library size) normalized
#'     counts.\cr
#'     See [DESeq2::sizeFactors()] for details.
#'   - `"fpkm"`: **F**ragments **p**er **k**ilobase per **m**illion mapped
#'     fragments.\cr
#'     Requires `fpkm = TRUE` in [bcbioRNASeq()] call and gene annotations in
#'     [`rowRanges()`][SummarizedExperiment::rowRanges] with defined
#'     [`width`][BiocGenerics::width].\cr
#'     See [DESeq2::fpkm()] for details.
#'   - `"vst"`: **V**ariance-**s**tabilizing **t**ransformation (log2).\cr
#'     Requires `vst = TRUE` to be set during [bcbioRNASeq()] call.\cr
#'     See `DESeq2::varianceStabilizingTransformation` for more information.
#'   - `"rlog"`: **R**egularized **log** transformation (log2).\cr
#'     Requires `rlog = TRUE` to be set during [bcbioRNASeq()] call.\cr
#'     See [DESeq2::rlog()] for details.
#'   - `"tmm"`: **T**rimmed **m**ean of **M**-values.\cr
#'     Calculated on the fly.\cr
#'     See [edgeR::calcNormFactors()] for details.
#'   - `"rle"`: **R**elative **l**og **e**xpression transformation.\cr
#'     Calculated on the fly.\cr
#'     See [relativeLogExpression()] for details.
#'
#' [featureCounts]: http://bioinf.wehi.edu.au/featureCounts/
#' [tximport]: https://bioconductor.org/packages/tximport/
#'
#' @return `matrix`.
#'
#' @references
#' - TMM: Robinson and Oshlack (2010).
#' - RLE: Anders and Huber (2010).
#'
#' @seealso
#' - `tximport::tximport()`.
#' - `DESeq2::counts()`.
#' - `DESeq2::sizeFactors()`.
#' - `DESeq2::varianceStabilizingTransformation()`.
#' - `DESeq2::rlog()`.
#' - `DESeq2::fpkm()`.
#' - `edgeR::calcNormFactors()`.
#' - `tmm()`.
#' - `relativeLogExpression()`.
#'
#' @examples
#' data(bcb)
#' summary(counts(bcb))
NULL



#' @rdname counts
#' @name counts
#' @importFrom BiocGenerics counts
#' @export
NULL

#' @rdname counts
#' @name counts<-
#' @importFrom BiocGenerics counts<-
#' @usage NULL
#' @export
NULL



## Updated 2019-07-23.
`counts,bcbioRNASeq` <-  # nolint
    function(object, normalized = FALSE) {
        validObject(object)
        assert(
            isAny(normalized, classes = c("character", "logical")),
            ## Ensure that primary assay matches counts.
            identical(assayNames(object)[[1L]], "counts")
        )
        ## Restrict the `normalized` arguments for transcript-level objects.
        if (.isTranscriptLevel(object)) {
            assert(isSubset(normalized, list(FALSE, "tpm")))
        }
        if (identical(normalized, FALSE)) {
            assayName <- "counts"
        } else if (
            identical(normalized, TRUE) ||
            ## `sf` is short for library size factor normalized.
            identical(normalized, "sf")
        ) {
            assayName <- "normalized"
        } else {
            assayName <- match.arg(arg = normalized, choices = normalizedCounts)
        }
        if (assayName == "tmm") {
            counts <- assays(object)[["counts"]]
            counts <- tmm(counts)
        } else if (assayName == "rle") {
            counts <- assays(object)[["counts"]]
            counts <- relativeLogExpression(counts)
        } else {
            if (!isSubset(assayName, assayNames(object))) {
                stop(sprintf("%s counts are not defined in object.", assayName))
            }
            counts <- assays(object)[[assayName]]
        }
        assert(is.matrix(counts))
        counts
    }



#' @rdname counts
#' @export
setMethod(
    f = "counts",
    signature = signature("bcbioRNASeq"),
    definition = `counts,bcbioRNASeq`
)
