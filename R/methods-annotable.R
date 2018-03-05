#' Ensembl Annotations
#'
#' @rdname annotable
#' @name annotable
#' @author Michael Steinbaugh
#'
#' @importFrom basejump annotable
#'
#' @inheritParams general
#'
#' @return [data.frame]
#'
#' @examples
#' load(system.file("extdata/bcb.rda", package = "bcbioRNASeq"))
#'
#' # bcbioRNASeq
#' annotable(bcb) %>% glimpse()
NULL



# Methods ======================================================================
#' @rdname annotable
#' @export
setMethod(
    "annotable",
    signature("bcbioRNASeq"),
    function(object) {
        .Deprecated("rowData")
        rowData(object)
    })
