#' `tximport` Accessor Shortcut
#'
#' @rdname txi
#' @name txi
#' @author Michael Steinbaugh
#'
#' @return [tximport::tximport()] list.
#'
#' @examples
#' txi(bcb) %>% names
NULL



# Methods ====
#' @rdname txi
#' @export
setMethod("txi", "bcbioRNADataSet", function(object) {
    bcbio(object, type = "tximport")
})
