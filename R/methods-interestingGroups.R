#' Interesting Groups
#'
#' @rdname interestingGroups
#' @name interestingGroups
#'
#' @importFrom bcbioBase interestingGroups interestingGroups<-
#'
#' @inheritParams AllGenerics
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
#' @importFrom bcbioBase checkInterestingGroups
#' @importFrom S4Vectors metadata
#' @export
setMethod(
    "interestingGroups<-",
    signature(object = "bcbioRNASeq", value = "character"),
    function(object, value) {
        sampleMetadata <- sampleMetadata(object)
        interestingGroups <- checkInterestingGroups(
            object = sampleMetadata,
            interestingGroups = value)
        metadata(object)[["interestingGroups"]] <- interestingGroups
        validObject(object)
        object
    })



#' @rdname interestingGroups
#' @importFrom S4Vectors metadata
#' @export
setMethod(
    "interestingGroups<-",
    signature(object = "bcbioRNASeq", value = "NULL"),
    function(object, value) {
        warning(paste(
            "'interestingGroups' is 'NULL'.",
            "Defaulting to 'sampleName'."
        ), call. = FALSE)
        metadata(object)[["interestingGroups"]] <- "sampleName"
        validObject(object)
        object
    })
#'
