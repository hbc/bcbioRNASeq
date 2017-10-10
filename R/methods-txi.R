#' `tximport` Accessor Shortcut
#'
#' @rdname txi
#' @name txi
#' @author Michael Steinbaugh
#' @keywords internal
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
setMethod(
    "txi",
    signature("bcbioRNASeqANY"),
    function(object) {
        bcbio(object, type = "tximport")
    })
