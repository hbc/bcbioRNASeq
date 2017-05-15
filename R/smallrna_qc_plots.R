#' Small RNA-seq quality control plots.
#'
#' @rdname smallrna_qc_plots
#' @author Lorena Pantano
#'
#' @param run bcbio-nextgen run.
#'
#' @return ggplot figure.



#' @rdname smallrna_qc_plots
#' @description Plot size distribution of small RNA-seq data
#' @export
plot_size_distribution <- function(run) {
    import_tidy_verbs()
    fns <- file.path(run$sample_dirs,
                     paste(names(run$sample_dirs),
                           "ready.trimming_stats",
                           sep = "-"))
    names(fns) <- names(run$sample_dirs)
    tab <- data.frame()
    for (sample in rownames(run$metadata)) {
        d <- read.table(fns[sample], sep = "")
        tab <- rbind(
            tab, d %>%
                mutate(sample = .data$sample,
                       group = run$metadata[.data$sample, run$intgroup]))
    }

    reads_adapter <- tab %>%
        group_by(!!!syms(c("sample", "group"))) %>%
        summarise(total = sum(.data$V2))

    ggdraw() +
        draw_plot(
            ggplot(reads_adapter,
                   aes_(x = ~sample, y = ~total, fill = ~group)) +
                geom_bar(stat = "identity", position = "dodge") +
                ggtitle("total number of reads with adapter") +
                ylab("# reads") +
                theme(
                    axis.text.x = element_text(
                        angle = 90, vjust = 0.5, hjust = 1)), 0, 0.5, 1, 0.5) +
        draw_plot(
            ggplot(tab, aes_(x = ~V1, y = ~V2, group = ~sample)) +
                geom_bar(stat = "identity", position = "dodge") +
                facet_wrap(~group, ncol = 2) +
                ggtitle("size distribution") +
                ylab("# reads") + xlab("size") +
                theme(
                    axis.text.x = element_text(
                        angle = 90, vjust = 0.5, hjust = 1)), 0, 0, 1, 0.5)
}



#' @rdname smallrna_qc_plots
#' @description Plot of total miRNA counts
#' @export
plot_mirna_counts <- function(run) {
    data.frame(sample = colnames(counts(run$srna_counts)),
               total = colSums(counts(run$srna_counts))) %>%
    ggplot +
        geom_bar(aes_(x = ~sample, y = ~total), stat = "identity") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
}
