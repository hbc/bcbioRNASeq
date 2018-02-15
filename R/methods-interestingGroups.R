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
#' load(system.file(
#'     file.path("extdata", "bcb.rda"),
#'     package = "bcbioRNASeq"))
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
#' @importFrom S4Vectors metadata
#' @export
setMethod(
    "interestingGroups",
    signature("bcbioRNASeq"),
    function(object) {
        metadata(object)[["interestingGroups"]]
    })



#' @rdname interestingGroups
#' @importFrom bcbioBase assert_formal_interesting_groups
#' @importFrom S4Vectors metadata
#' @export
setMethod(
    "interestingGroups<-",
    signature(object = "bcbioRNASeq", value = "character"),
    function(object, value) {
        assert_formal_interesting_groups(
            object = sampleMetadata(object),
            interestingGroups = value)
        metadata(object)[["interestingGroups"]] <- value
        validObject(object)
        object
    })
