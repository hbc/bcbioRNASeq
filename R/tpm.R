#' Transcripts per million
#'
#' Save TPM values from a tximport counts object
#'
#' Links:
#'
#' - [tximport vignette](https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html)
#' - [salmon getting started guide](https://combine-lab.github.io/salmon/getting_started/)
#' - [Bioconductor thread](https://support.bioconductor.org/p/84883/#84929)
#'
#' @author Michael Steinbaugh
#'
#' @param txi tximport list object
#'
#' @return TPM (transcripts per million) matrix
#' @export
tpm <- function(txi) {
    if (!identical(
        names(txi),
        c("abundance",
          "counts",
          "length",
          "countsFromAbundance")
    )) {
        stop("tximport list is required")
    }

    tpm <- txi$abundance

    return(tpm)
}
