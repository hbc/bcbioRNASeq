#' Plot Dispersion Estimates
#'
#' [DESeq2::plotDispEsts()] wrapper supporting a [bcbioRNADataSet].
#'
#' @rdname plotDispersion
#' @name plotDispersion
#'
#' @seealso [DESeq2::plotDispEsts()].
#'
#' @return Dispersion plot [ggplot].
#' @export
#'
#' @examples
#' data(bcb, dds)
#'
#' # bcbioRNADataSet
#' plotDispersion(bcb)
#'
#' # DESeqDataSet
#' plotDispersion(dds)
NULL



# Constructors ====
.plotDispersion <- function(object) {
    DESeq2::plotDispEsts(object)
}



# Methods ====
#' @rdname plotDispersion
#' @export
setMethod("plotDispersion", "bcbioRNADataSet", function(object) {
    bcbio(object, "DESeqDataSet") %>% .plotDispersion
})



#' @rdname plotDispersion
#' @export
setMethod("plotDispersion", "DESeqDataSet", .plotDispersion)
