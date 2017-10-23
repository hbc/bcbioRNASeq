#' Interesting Groups
#'
#' @rdname interestingGroups
#' @name interestingGroups
#'
#' @importFrom basejump interestingGroups
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
