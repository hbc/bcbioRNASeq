#' Interesting Groups
#'
#' @rdname interestingGroups
#' @name interestingGroups
#'
#' @inheritParams AllGenerics
#'
#' @return Character vector.
#'
#' @examples
#' data(bcb)
#' interestingGroups(bcb)
NULL



# Constructors ====
.interestingGroups <- function(object) {
    metadata(object)[["interestingGroups"]]
}



# Methods ====
#' @rdname interestingGroups
#' @export
setMethod(
    "interestingGroups",
    signature = "bcbioRNASeqANY",
    definition = .interestingGroups)
