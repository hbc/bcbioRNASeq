#' Slot aligned counts into object
#'
#' This function loads aligned counts (e.g. STAR, HISAT2) from the bcbio final
#' output directory into the bcbioRNASeq object, so we can visually inspect
#' correlations with the primary pseudoaligned counts.
#'
#' @name slotAlignedCounts
#' @note Updated 2019-08-07.
#'
#' @inheritParams acidroxygen::params
#'
#' @examples
#' uploadDir <- system.file("extdata/bcbio", package = "bcbioRNASeq")
#' ## Fast mode skips import of aligned counts.
#' bcb <- bcbioRNASeq(uploadDir, fast = TRUE)
#' "aligned" %in% assayNames(bcb)
#' bcb <- slotAlignedCounts(bcb)
#' "aligned" %in% assayNames(bcb)
NULL



## Updated 2019-07-23.
`slotAlignedCounts,bcbioRNASeq` <-  # nolint
    function(object) {
        validObject(object)
        assert(areDisjointSets(assayNames(object), "aligned"))
        assays(object)[["aligned"]] <- .featureCounts(
            projectDir = metadata(object)[["projectDir"]],
            samples = colnames(object),
            genes = rownames(object)
        )
        validObject(object)
        object
    }



#' @rdname slotAlignedCounts
#' @export
setMethod(
    f = "slotAlignedCounts",
    signature = signature("bcbioRNASeq"),
    definition = `slotAlignedCounts,bcbioRNASeq`
)
