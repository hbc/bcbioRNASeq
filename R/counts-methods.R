#' @name counts
#' @author Michael Steinbaugh, Lorena Pantano
#' @inherit AcidGenerics::counts
#' @note Updated 2021-02-22.
#'
#' @details
#' By default, return the raw counts. Normalized counts in a variety of formats
#' can be accessed using the `normalized` argument.
#'
#' @inheritParams AcidRoxygen::params
#' @param normalized `character(1)` or `logical(1)`.
#'   Normalization method to apply:
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
#'     See `DESeq2::sizeFactors` for details.
#'   - `"fpkm"`: **F**ragments **p**er **k**ilobase per **m**illion mapped
#'     fragments.\cr
#'     Requires `fast = FALSE` in [bcbioRNASeq()] call and gene annotations in
#'     [`rowRanges()`][SummarizedExperiment::rowRanges] with defined
#'     [`width`][BiocGenerics::width].\cr
#'     See [DESeq2::fpkm()] for details.
#'   - `"vst"`: **V**ariance-**s**tabilizing **t**ransformation (log2).\cr
#'     Requires `fast = FALSE` to be set during [bcbioRNASeq()] call.\cr
#'     See `DESeq2::varianceStabilizingTransformation` for more information.
#'   - `"tmm"`: **T**rimmed **m**ean of **M**-values.\cr
#'     Calculated on the fly.\cr
#'     See [edgeR::calcNormFactors()] for details.
#'   - `"rle"`: **R**elative **l**og **e**xpression transformation.\cr
#'     Calculated on the fly.\cr
#'     See [relativeLogExpression()] for details.
#'   - `"rlog"`: *Deprecated*.
#'     **R**egularized **log** transformation (log2).\cr
#'     No longer calculated automatically during [bcbioRNASeq()] call, but may
#'     be defined in legacy objects.\cr
#'     See [DESeq2::rlog()] for details.\cr
#'     Note that VST is more performant and now recommended by default instead.
#'
#'     Note that `logical(1)` support only applies to `counts()`. Other
#'     functions in the package require `character(1)` and use
#'     [`match.arg()`][base::match.arg] internally.
#'
#' [featureCounts]: http://bioinf.wehi.edu.au/featureCounts/
#' [tximport]: https://bioconductor.org/packages/tximport/
#'
#' @param value Value to assign.
#' @param ... Additional arguments.
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



## Updated 2020-09-23.
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
            isTRUE(normalized) ||
            ## `sf` is short for library size factor normalized.
            identical(normalized, "sf")
        ) {
            assayName <- "normalized"
        } else {
            assayName <- match.arg(arg = normalized, choices = .normalized)
        }
        if (identical(assayName, "tmm")) {
            counts <- assays(object)[["counts"]]
            counts <- tmm(counts)
        } else if (identical(assayName, "rle")) {
            counts <- assays(object)[["counts"]]
            counts <- relativeLogExpression(counts)
        } else {
            if (!isSubset(assayName, assayNames(object))) {
                abort(sprintf(
                    "{.var %s} counts are not defined in object.", assayName
                ))
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
