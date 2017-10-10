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



# Methods ====
#' @rdname interestingGroups
#' @export
setMethod(
    "interestingGroups",
    signature("bcbioRNASeqANY"),
    function(object) {
        metadata(object)[["interestingGroups"]]
    })
