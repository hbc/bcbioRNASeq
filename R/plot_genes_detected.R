#' Genes detected plot
#'
#' @author Michael Steinbaugh
#'
#' @import ggplot2
#'
#' @param summary \code{bcbio-rnaseq} summary report
#' @param counts Raw counts matrix
#'
#' @export
#'
#' @examples
#' \dontrun{
#' plot_genes_detected(summary, raw_counts)
#' }
plot_genes_detected <- function(summary, counts) {
    summary %>%
        ggplot2::ggplot(
            ggplot2::aes_(x = ~description,
                          y = colSums(counts > 0),
                          fill = ~intgroup)
        ) +
        ggplot2::ggtitle("genes detected") +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::geom_hline(color = "green",
                            size = 2,
                            yintercept = 20000) +
        ggplot2::xlab("sample") +
        ggplot2::ylab("gene count") +
        ggplot2::coord_flip()
}
