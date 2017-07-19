#' `tximport` Accessor Shortcut
#'
#' @rdname txi
#' @author Michael Steinbaugh
#'
#' @return [tximport::tximport()] list.
#' @export
setMethod("txi", "bcbioRNADataSet", function(object) {
    bcbio(object, type = "tximport")
})
