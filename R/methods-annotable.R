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
#' annotable(bcb) %>% str()
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
