#' Small RNA-seq quality control plots
#'
#' @rdname plots-smallrna
#' @author Lorena Pantano
#' @author Michael Steinbaugh
#'
#' @param run [bcbioRnaDataSet].
#' @return [ggplot].



#' @rdname plots-smallrna
#' @description Plot size distribution of small RNA-seq data.
#' @export
plot_size_distribution <- function(run) {
    run <- metadata(run)
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
                mutate(sample = sample,
                       group = run$metadata[sample, run$interesting_groups]))
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
                        angle = 90, vjust = 0.5, hjust = 1)),
            0, 0.5, 1, 0.5) +
        draw_plot(
            ggplot(tab, aes_(x = ~V1, y = ~V2, group = ~sample)) +
                geom_bar(stat = "identity", position = "dodge") +
                facet_wrap(~group, ncol = 2) +
                ggtitle("size distribution") +
                ylab("# reads") + xlab("size") +
                theme(
                    axis.text.x = element_text(
                        angle = 90, vjust = 0.5, hjust = 1)),
            0, 0, 1, 0.5)
}



#' @rdname plots-smallrna
#' @description Plot of total miRNA counts.
#' @export
plot_mirna_counts <- function(run) {
    # [TODO] `txi()` accessor doesn't exist
    .t <- data.frame(sample = colnames(txi(run)),
                     total = colSums(txi(run)))
    cs <- apply(txi(run), 2L, function(x) {
        cumsum(sort(x, decreasing = TRUE))
    }) %>% as.data.frame
    cs$pos <- 1:nrow(cs)
    plot_grid(
        ggplot(.t) +
            geom_bar(aes_(x = ~sample,
                          y = ~total),
                     stat = "identity") +
            theme(axis.text.x = element_text(
                angle = 90L, hjust = 1L, vjust = 0.5)),
        ggplot(melt(txi(run))) +
            geom_boxplot(aes_(x = ~X2, y = ~value)) +
            xlab("") +
            ylab("expression") +
            scale_y_log10() +
            theme(axis.text.x = element_text(
                angle = 90L, hjust = 1L, vjust = 0.5)),
        ggplot((melt(cs, id.vars = "pos"))) +
            geom_line(aes_(x = ~pos,
                           y = ~value,
                           color = ~variable)) +
            xlim(0L, 50L) +
            scale_y_log10(),
        nrow = 3L)
}


#' @rdname plots-smallrna
#' @description Clustering small RNA samples.
#' @export
plot_srna_clusters <- function(run) {
    counts <- bcbio(run)
    design <- metadata(run)$metadata
    dds <- DESeqDataSetFromMatrix(
        counts[rowSums(counts > 0L) > 3L, ],
        colData = design,
        design = ~1)
    vst <- rlog(dds, betaPriorVar = FALSE)
    annotation_col <- design %>%
        .[, metadata(run)[["interesting_groups"]], drop = FALSE]
    pheatmap(assay(vst),
             annotation_col = annotation_col,
             clustering_distance_cols = "correlation",
             clustering_method = "ward.D",
             scale = "row",
             show_rownames = FALSE)
    plot_pca(metadata(run), vst)
}
