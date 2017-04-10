#' Transcripts per million
#'
#' Save TPM values from a tximport counts object
#'
#' @author Michael Steinbaugh
#'
#' @param txi \code{tximport} object
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
