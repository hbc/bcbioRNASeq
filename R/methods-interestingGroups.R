#' Interesting Groups
#'
#' @rdname interestingGroups
#' @name interestingGroups
#'
#' @importFrom bcbioBase interestingGroups interestingGroups<-
#'
#' @inheritParams general
#'
#' @return Character vector.
#'
#' @examples
#' load(system.file("extdata/bcb.rda", package = "bcbioRNASeq"))
#'
#' # bcbioRNASeq
#' interestingGroups(bcb)
#'
#' # Assignment support
#' interestingGroups(bcb) <- "sampleID"
#' interestingGroups(bcb)
NULL



# Methods ======================================================================
#' @rdname interestingGroups
#' @export
setMethod(
    "interestingGroups",
    signature("bcbioRNASeq"),
    function(object) {
        metadata(object)[["interestingGroups"]]
    })



#' @rdname interestingGroups
#' @export
setMethod(
    "interestingGroups<-",
    signature(object = "bcbioRNASeq", value = "character"),
    function(object, value) {
        assertFormalInterestingGroups(
            sampleMetadata(object),
            interestingGroups = value)
        metadata(object)[["interestingGroups"]] <- value
        validObject(object)
        object
    })
