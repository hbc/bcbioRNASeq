#' Interesting Groups
#'
#' @rdname interestingGroups
#' @name interestingGroups
#'
#' @inheritParams AllGenerics
#'
#' @return Character vector.
#' @export
setMethod(
    "interestingGroups",
    signature = "bcbioRNASeqANY",
    definition = function(object) {
        metadata(object)[["interestingGroups"]]
    })
