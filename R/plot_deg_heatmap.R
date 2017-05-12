#' Differentially expressed gene (DEG) heatmap
#'
#' @author Michael Steinbaugh
#'
#' @param run bcbio-nextgen run.
#' @param dds \linkS4class{DESeqDataSet}.
#' @param contrast Contrast.
#' @param alpha Alpha level cutoff.
#' @param lfc Log fold change ratio (base 2) cutoff.
#'
#' @return Gene heatmap.
#' @export
plot_deg_heatmap <- function(
    run,
    dds,
    contrast,
    alpha = 0.1,
    lfc = 0) {
    check_run(run)
    check_dds(dds)

    # Annotation metadata
    annotation <- colData(dds) %>% as.data.frame %>% .[, run$intgroup]

    # DE genes, used for subsetting the rlog counts matrix
    res <- results(dds,
                   contrast = contrast,
                   alpha = alpha)

    # Output results to data frame and subset by alpha and lfc cutoffs
    res_df <- res %>%
        as.data.frame %>%
        rownames_to_column("ensembl_gene_id") %>%
        as_tibble %>%
        .[.$padj < alpha, ] %>%
        .[which(.$log2FoldChange > lfc | .$log2FoldChange < -lfc), ]

    # rlog transform and subset the counts
    counts <- dds %>% rlog %>% assay %>% .[res_df$ensembl_gene_id, ]

    if (nrow(counts) <= 100) {
        # Change rownames to readable external gene names
        show_rownames <- TRUE
        rownames(counts) <- run$ensembl[rownames(counts), "external_gene_name"]
    } else {
        show_rownames <- FALSE
    }

    name <- deparse(substitute(dds))
    contrast_name <- res_contrast_name(res)

    pheatmap(counts,
             annotation = annotation,
             clustering_distance_cols = "correlation",
             clustering_method = "ward.D2",
             main = paste(name, contrast_name, sep = " : "),
             scale = "row",
             show_rownames = show_rownames)
}
