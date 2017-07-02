#' Access tximport data in bcbioRNADataSet
#'
#' @rdname txi
#' @author Michael Steinbaugh
#'
#' @param object Object.
#'
#' @return [tximport::tximport()] list.
#' @export
setMethod("txi", "bcbioRNADataSet", function(object) {
    bcbio(object, type = "tximport")
})
