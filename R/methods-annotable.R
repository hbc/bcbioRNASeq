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
#' load(system.file(
#'     file.path("extdata", "bcb.rda"),
#'     package = "bcbioRNASeq"))
#'
#' # bcbioRNASeq
#' annotable(bcb) %>% glimpse()
NULL



# Constructors =================================================================
.annotable.bcbioRNASeq <- function(object) {  # nolint
    data <- rowData(object)
    rownames(data) <- slot(object, "NAMES")
    as.data.frame(data)
}



# Methods ======================================================================
#' @rdname annotable
#' @importFrom S4Vectors metadata
#' @export
setMethod(
    "annotable",
    signature("bcbioRNASeq"),
    .annotable.bcbioRNASeq)
