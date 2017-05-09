#' Small RNA-seq quality control plots
#'
#' @rdname smallrna_qc_plots
#' @author Lorena Pantano
#'
#' @param run \code{bcbio-nextgen} run
#'
#' @return ggplot figure



#' @rdname smallrna_qc_plots
#' @description Plot size distribution of small RNA-seq data
#' @export
plot_size_distribution <- function(run) {
    # [fix] NSE handling for group, total, V1, V2 using dplyr 0.6.0
    fns <- file.path(run$sample_dirs,
                     paste(names(run$sample_dirs),
                           "ready.trimming_stats",
                           sep = "-"))
    names(fns) <- names(run$sample_dirs)

    tab <- data.frame()
    for (sample in rownames(run$metadata)) {
        d <- read.table(fns[sample], sep = "")
        # [fix] `mutate` NSE - use `.data$var`
        tab <- rbind(
            tab, d %>% mutate(sample = sample,
                              group = run$metadata[sample, run$intgroup]))
    }

    reads_adapter <- tab %>%
        # [fix] `group_by`, `summarize` NSE -- use `.data$var` or `!!sym("var")`
        group_by(sample, group) %>%
        summarize(total = sum(V2))

    # [fix] `aes` to `aes_` NSE -- use `~var`
    ggdraw() +
        draw_plot(
            ggplot(reads_adapter, aes(x = sample,y = total, fill = group)) +
                geom_bar(stat = "identity", position = "dodge") +
                ggtitle("total number of reads with adapter") +
                ylab("# reads") +
                theme(
                    axis.text.x = element_text(angle = 90,
                                               vjust = 0.5,
                                               hjust = 1)), 0, 0.5, 1, 0.5) +
        draw_plot(
            ggplot(tab, aes(x = V1, y = V2, group = sample)) +
                geom_bar(stat = "identity", position = "dodge") +
                facet_wrap(~group, ncol = 2) +
                ggtitle("size distribution") +
                ylab("# reads") + xlab("size") +
                theme(
                    axis.text.x = element_text(angle = 90,
                                               vjust = 0.5,
                                               hjust = 1)), 0, 0, 1, 0.5)
}

#' @rdname smallrna_qc_plots
#' @description Plot of total miRNA counts
#' @export
plot_mirna_counts <- function(run){
    ggplot( data.frame(sample=colnames(counts(run$srna_counts)),
                       total=colSums(counts(run$srna_counts)))) +
        geom_bar(aes(x=sample,y=total), stat='identity')+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}
