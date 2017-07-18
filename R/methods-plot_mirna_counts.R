#' Plot of total miRNA counts
#'
#' @rdname plot_mirna_counts
#' @author Lorena Pantano
#'
#' @return [ggplot].
#' @export
setMethod("plot_mirna_counts", "bcbioRNADataSet", function(object) {
    counts <- counts(object)
    .t <- data.frame(sample = colnames(counts),
                     total = colSums(counts))
    cs <- apply(counts, 2L, function(x) {
        cumsum(sort(x, decreasing = TRUE))
    }) %>% as.data.frame
    cs[["pos"]] <- 1L:nrow(cs)
    plot_grid(
        ggplot(.t) +
            geom_bar(aes_(x = ~sample,
                          y = ~total),
                     stat = "identity") +
            theme(axis.text.x = element_text(
                angle = 90L, hjust = 1L, vjust = 0.5)),
        ggplot(melt(counts)) +
            geom_boxplot(aes_(x = ~X2, y = ~value)) +
            xlab("") +
            ylab("expression") +
            scale_y_log10() +
            theme(axis.text.x = element_text(
                angle = 90L, hjust = 1L, vjust = 0.5)),
        ggplot(melt(cs, id.vars = "pos")) +
            geom_line(aes_(x = ~pos,
                           y = ~value,
                           color = ~variable)) +
            xlim(0L, 50L) +
            scale_y_log10(),
        nrow = 3L)
})
