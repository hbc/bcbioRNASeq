#' Differentially expressed gene heatmap
#'
#' @author Michael Steinbaugh
#'
#' @param run bcbio-nextgen run.
#' @param res \linkS4class{DESeqResults}.
#' @param dt \linkS4class{DESeqTransform}. We recommend passing in [rlog()]
#'   counts in as the transform object by default.
#' @param lfc Log2 fold change ratio cutoff.
#'
#' @return Gene heatmap.
#' @export
plot_deg_heatmap <- function(
    run,
    res,
    dt,
    lfc = 0) {
    check_run(run)
    check_res(res)
    check_dt(dt)
    import_tidy_verbs()

    alpha <- res@metadata$alpha
    metadata <- intgroup_as_factor(run)

    # Output results to data frame and subset by alpha and lfc cutoffs
    res_df <- res %>%
        as.data.frame %>%
        rownames_to_column("ensembl_gene_id") %>%
        as_tibble %>%
        filter(.data$padj < alpha) %>%
        .[which(.$log2FoldChange > lfc | .$log2FoldChange < -lfc), ]

    # rlog transform and subset the counts
    counts <- dt %>% assay %>% .[res_df$ensembl_gene_id, ]

    if (nrow(counts) <= 100) {
        # Change rownames to readable external gene names
        show_rownames <- TRUE
        gene_names <- gene_level_annotations(run) %>%
            select(!!!syms(c("ensembl_gene_id", "external_gene_name"))) %>%
            as.data.frame %>%
            set_rownames(.$ensembl_gene_id)
        rownames(counts) <- gene_names[rownames(counts), "external_gene_name"]
    } else {
        show_rownames <- FALSE
    }

    res_name <- deparse(substitute(res))
    contrast_name <- res_contrast_name(res)
    dt_name <- deparse(substitute(dt))

    pheatmap(counts,
             annotation = metadata,
             clustering_distance_cols = "correlation",
             clustering_method = "ward.D2",
             main = paste(res_name, contrast_name, dt_name, sep = " : "),
             scale = "row",
             show_rownames = show_rownames)
}
