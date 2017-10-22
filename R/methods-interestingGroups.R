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
#' data(bcb)
#' interestingGroups(bcb)
NULL



# Methods ====
#' @rdname interestingGroups
#' @export
setMethod(
    "interestingGroups",
    signature("bcbioRNASeq"),
    function(object) {
        metadata(object)[["interestingGroups"]]
    })
