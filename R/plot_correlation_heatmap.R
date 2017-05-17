#' Correlation matrix heatmap
#'
#' @author Michael Steinbaugh
#'
#' @param run bcbio-nextgen run.
#' @param dt \linkS4class{DESeqTransform} generated from [DESeq2::rlog()] or
#'   [DESeq2::vst()] on a \linkS4class{DESeqDataSet}.
#' @param method Correlation coefficient (or covariance) to be computed.
#'   Defaults to `pearson` but `spearman` can also be used.
#'
#' @return Correlation heatmap.
#' @export
plot_correlation_heatmap <- function(
    run,
    dt,
    method = "pearson") {
    check_run(run)
    check_dt(dt)
    if (!method %in% c("pearson", "spearman")) {
        stop("Support methods: pearson, spearman")
    }
    name <- deparse(substitute(dt))

    # Get counts and annotations from DESeqTransform object
    counts <- assay(dt)
    annotation <- colData(dt) %>%
        .[, run$intgroup] %>%
        as.data.frame

    # Pearson or Spearman correlation methods are supported
    if (!method %in% c("pearson", "spearman")) {
        stop("invalid correlation regression method.
             must use pearson or spearman.")
    }

    counts %>%
        cor(method = method) %>%
        pheatmap(
            annotation = annotation,
            main = paste(
                "correlation",
                method,
                name,
                sep = " : "),
            show_colnames = TRUE,
            show_rownames = TRUE)
}
