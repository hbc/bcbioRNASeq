#' `tximport` Accessor Shortcut
#'
#' @rdname txi
#' @name txi
#' @author Michael Steinbaugh
#'
#' @inheritParams AllGenerics
#'
#' @return [tximport::tximport()] list.
#'
#' @examples
#' txi(bcb) %>%
#'     names()
NULL



# Methods ====
#' @rdname txi
#' @export
setMethod("txi", "bcbioRNASeqANY", function(object) {
    bcbio(object, type = "tximport")
})
