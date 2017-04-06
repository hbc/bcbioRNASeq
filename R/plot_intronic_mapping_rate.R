#' Intronic mapping rate plot
#'
#' @author Michael Steinbaugh
#'
#' @import ggplot2
#'
#' @param bcbio bcbio list object
#'
#' @export
#'
#' @examples
#' \dontrun{
#' plot_intronic_mapping_rate(bcbio)
#' }
plot_intronic_mapping_rate <- function(bcbio) {
    check_bcbio(bcbio)
    plot <- import_summary(bcbio) %>%
        ggplot2::ggplot(
            ggplot2::aes_(x = ~description,
                          # Multiple by 100 here for percentage
                          y = ~intronic_rate * 100,
                          fill = ~qc_color)
        ) +
        ggplot2::ggtitle("intronic mapping rate") +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::geom_hline(color = "orange",
                            size = 2,
                            yintercept = 20) +
        ggplot2::labs(x = "sample",
                      y = "intronic mapping rate (%)",
                      fill = "") +
        ggplot2::ylim(0, 100) +
        ggplot2::coord_flip()
    print(plot)
}
