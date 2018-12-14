#' @name params
#' @inherit basejump::params
#' @keywords internal
#' @param normalized `character(1)`.
#'   Which normalization method to apply:
#'
#'   - "`tpm`": Transcripts per million (tximport).
#'   - "`tmm`": edgeR trimmed mean of M-values. Calculated on the fly.
#'   - "`rlog`": DESeq2 **log2** regularized log transformation.
#'   - "`vst`": DESeq2 **log2** variance stabilizing transformation.
NULL
