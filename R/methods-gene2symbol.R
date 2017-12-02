#' Gene to Symbol Mappings
#'
#' @rdname gene2symbol
#' @name gene2symbol
#' @author Michael Steinbaugh
#'
#' @importFrom basejump gene2symbol
#'
#' @inheritParams AllGenerics
#'
#' @examples
#' bcb <- examples[["bcb"]]
#' gene2symbol(bcb) %>% head()
NULL



# Constructors ====
.gene2symbol <- function(object) {
    annotable <- annotable(bcb)
    if (is.null(annotable)) {
        return(NULL)
    }
    cols <- c("ensgene", "symbol")
    if (!all(cols %in% colnames(annotable))) {
        stop(paste(
            toString(cols),
            "missing from internal annotable"
        ), call. = FALSE)
    }
    annotable[, cols]
}


# Methods ====
#' @rdname gene2symbol
#' @export
setMethod(
    "gene2symbol",
    signature("bcbioRNASeq"),
    .gene2symbol)
