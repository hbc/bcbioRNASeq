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
#'
#' @examples
#' bcb <- examples[["bcb"]]
#' annotable(bcb) %>%
#'     glimpse()
NULL



# Methods ====
#' @rdname annotable
#' @importFrom S4Vectors metadata
#' @export
setMethod(
    "annotable",
    signature("bcbioRNASeq"),
    function(object) {
        annotable <- as.data.frame(rowData(object))
        rownames(annotable) <- slot(object, "NAMES")
        annotable
    })
