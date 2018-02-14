#' Gene to Symbol Mappings
#'
#' @rdname gene2symbol
#' @name gene2symbol
#' @author Michael Steinbaugh
#'
#' @importFrom basejump gene2symbol
#'
#' @inheritParams general
#'
#' @return [data.frame] containing Ensembl gene identifier (`ensgene`) and
#'   symbol (`symbol`) mappings.
#'
#' @examples
#' load(system.file(
#'     file.path("extdata", "bcb.rda"),
#'     package = "bcbioRNASeq"))
#'
#' # bcbioRNASeq
#' gene2symbol(bcb) %>% head()
NULL



# Constructors =================================================================
.gene2symbol.bcbioRNASeq <- function(object) {  # nolint
    annotable <- annotable(object)
    assert_is_data.frame(annotable, severity = "warning")
    if (is.null(annotable)) {
        return(invisible())
    }
    assert_is_subset(c("ensgene", "symbol"), colnames(annotable))
    annotable[, cols, drop = FALSE]
}


# Methods ======================================================================
#' @rdname gene2symbol
#' @export
setMethod(
    "gene2symbol",
    signature("bcbioRNASeq"),
    .gene2symbol.bcbioRNASeq)
