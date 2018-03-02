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
#' load(system.file("extdata/bcb.rda", package = "bcbioRNASeq"))
#'
#' # bcbioRNASeq
#' gene2symbol(bcb) %>% head()
NULL



# Methods ======================================================================
#' @rdname gene2symbol
#' @export
setMethod(
    "gene2symbol",
    signature("bcbioRNASeq"),
    function(object) {
        x <- rowData(object)
        cols <- c("ensgene", "symbol")
        assert_is_subset(cols, colnames(x))
        x <- x[, cols, drop = FALSE]
        assertIsGene2symbol(x)
        x
    })
