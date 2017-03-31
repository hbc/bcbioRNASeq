# https://support.bioconductor.org/p/84883/#84929


#' Transcripts per million
#'
#' @author Michael Steinbaugh
#'
#' @param txi tximport list object
#'
#' @return TPM (transcripts per million) data frame
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
