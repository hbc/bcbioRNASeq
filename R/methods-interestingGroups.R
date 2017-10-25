#' Interesting Groups
#'
#' @rdname interestingGroups
#' @name interestingGroups
#'
#' @importFrom basejump interestingGroups interestingGroups<-
#'
#' @inheritParams AllGenerics
#'
#' @return Character vector.
#'
#' @examples
#' interestingGroups(bcb)
NULL



# Methods ====
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
#' @importFrom S4Vectors metadata
#' @export
setMethod(
    "interestingGroups<-",
    signature(object = "bcbioRNASeq", value = "character"),
    function(object, value) {
        # Check the interesting groups against the sample metadata. `colData()`
        # can also be used here in place of `sampleMetadata()`, but
        # `sampleMetadata()` is used here for better code consistency with
        # the bcbioSingleCell package.
        sampleMetadata <- sampleMetadata(object)
        interestingGroups <- .checkInterestingGroups(
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
