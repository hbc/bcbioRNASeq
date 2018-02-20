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
        data <- rowData(object)
        assert_is_non_empty(data)
        rownames(data) <- slot(object, "NAMES")
        as.data.frame(data)
    })
