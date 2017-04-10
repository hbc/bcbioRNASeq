#' @rdname qc_plots
#' @description Genes detected plot
#'
#' @param raw_counts Raw counts matrix. Can be obtained from \code{DESeqDataSet}
#'   by running \code{counts(normalized = FALSE)}.
#'
#' @return Bar plot
#' @export
plot_genes_detected <- function(bcbio, raw_counts) {
    check_bcbio(bcbio)
    plot <- import_summary(bcbio) %>%
        ggplot(
            aes_(x = ~description,
                 y = colSums(raw_counts > 0),
                 fill = ~qc_color)
        ) +
        ggtitle("genes detected") +
        geom_bar(stat = "identity") +
        geom_hline(color = "green",
                   size = 2,
                   yintercept = 20000) +
        labs(x = "sample",
             y = "gene count",
             fill = "") +
        coord_flip()

    return(plot)
}
