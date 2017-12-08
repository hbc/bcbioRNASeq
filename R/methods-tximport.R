#' Access tximport Data
#'
#' @rdname tximport
#' @name tximport
#' @author Michael Steinbaugh
#' @keywords internal
#'
#' @inheritParams AllGenerics
#'
#' @return [tximport::tximport()] list.
#'
#' @examples
#' load(system.file(
#'     file.path("inst", "extdata", "bcb.rda"),
#'     package = "bcbioRNASeq"))
#'
#' # bcbioRNASeq
#' tximport(bcb) %>% names()
NULL



# Methods ====
#' @rdname tximport
#' @export
setMethod(
    "tximport",
    signature("bcbioRNASeq"),
    function(object) {
        bcbio(object, type = "tximport")
    })
