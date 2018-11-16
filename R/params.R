#' @name params
#' @inherit basejump::params
#' @keywords internal
#' @param lfcShrink `boolean`. Use shrunken log2 fold change (LFC) values.
#' @param normalized `string`. Which normalization method to apply:
#'   - "`tpm`": Transcripts per million (tximport).
#'   - "`tmm`": edgeR trimmed mean of M-values. Calculated on the fly.
#'   - "`rlog`": DESeq2 **log2** regularized log transformation.
#'   - "`vst`": DESeq2 **log2** variance stabilizing transformation.
#' @param results `scalar`. Position or name of `DESeqResults`.
NULL
