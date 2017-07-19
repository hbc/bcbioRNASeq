#' Plot Size Distribution of Small RNA-Seq Data
#'
#' @rdname plot_size_distribution
#' @author Lorena Pantano, Michael Steinbaugh
#' @family Small RNA-Seq Utilities
#'
#' @return [ggplot].
#' @export
setMethod("plot_size_distribution", "bcbioRNADataSet", function(object) {
    meta <- metadata(object)
    fns <- file.path(meta[["sample_dirs"]],
                     paste(names(meta[["sample_dirs"]]),
                           "ready.trimming_stats",
                           sep = "-"))
    names(fns) <- names(meta[["sample_dirs"]])
    tab <- data.frame()
    for (sample in rownames(meta[["metadata"]])) {
        d <- read.table(fns[sample], sep = "")
        tab <- rbind(
            tab, d %>%
                mutate(
                    sample = sample,
                    group = meta[["metadata"]][sample,
                                               meta[["interesting_groups"]]]))
    }

    reads_adapter <- tab %>%
        group_by(!!!syms(c("sample", "group"))) %>%
        summarise(total = sum(.data[["V2"]]))

    ggdraw() +
        draw_plot(
            ggplot(reads_adapter,
                   aes_(x = ~sample, y = ~total, fill = ~group)) +
                geom_bar(stat = "identity", position = "dodge") +
                ggtitle("total number of reads with adapter") +
                ylab("# reads") +
                theme(
                    axis.text.x = element_text(
                        angle = 90L, vjust = 0.5, hjust = 1L)),
            0L, 0.5, 1L, 0.5) +
        draw_plot(
            ggplot(tab, aes_(x = ~V1, y = ~V2, group = ~sample)) +
                geom_bar(stat = "identity", position = "dodge") +
                facet_wrap(~group, ncol = 2L) +
                ggtitle("size distribution") +
                ylab("# reads") + xlab("size") +
                theme(
                    axis.text.x = element_text(
                        angle = 90L, vjust = 0.5, hjust = 1L)),
            0L, 0L, 1L, 0.5)
})
