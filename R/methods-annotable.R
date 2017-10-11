#' Ensembl Annotations
#'
#' @rdname annotable
#' @name annotable
#' @author Michael Steinbaugh
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
    signature("bcbioRNASeqANY"),
    function(object) {
        metadata(object)[["annotable"]]
    })
