#' Ensembl Annotations
#'
#' @rdname annotable
#' @name annotable
#' @author Michael Steinbaugh
#'
#' @importFrom basejump annotable
#'
#' @inheritParams AllGenerics
#'
#' @return [data.frame]
NULL



# Methods ====
#' @rdname annotable
#' @export
setMethod(
    "annotable",
    signature("bcbioRNASeq"),
    function(object) {
        metadata(object)[["annotable"]]
    })
