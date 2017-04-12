#' Differentially expressed gene (DEG) heatmap
#'
#' @author Michael Steinbaugh
#'
#' @param run \code{bcbio-nextgen} run object
#' @param dds \code{DESeqDataSet}
#' @param contrast Contrast
#' @param alpha Alpha level cutoff
#' @param lfc Log fold change ratio (base 2) cutoff
#'
#' @return Gene heatmap
#' @export
deg_heatmap <- function(
    run,
    dds,
    contrast,
    alpha = 0.1,
    lfc = 0) {
    check_run(run)
    check_dds(dds)

    # rlog transform the counts
    rld <- rlog(dds)

    # Annotation metadata
    annotation <- dds %>%
        colData %>%
        as.data.frame %>%
        .[, run$intgroup]

    # DE genes, used for subsetting the rlog counts matrix
    res <- results(dds,
                   contrast = contrast,
                   alpha = alpha)

    # Output results to dataframe and subset by alpha and lfc cutoffs
    res_df <- res %>%
        as.data.frame %>%
        rownames_to_column("ensembl_gene_id") %>%
        filter_(.dots = ~padj < alpha) %>%
        filter_(.dots = ~log2FoldChange < -lfc | log2FoldChange > lfc)

    # Subset the transformed count matrix
    mat <- rld %>%
        assay %>%
        .[res_df$ensembl_gene_id, ]

    if (nrow(mat) <= 100) {
        show_rownames <- TRUE
        # Change rownames to readable external gene names
        gene_names <- ensembl_annotations(run, values = rownames(mat))
        if (identical(rownames(mat), rownames(gene_names))) {
            rownames(mat) <- gene_names$external_gene_name
        } else {
            stop("gene identifier rownames mismatch")
        }
    } else {
        show_rownames <- FALSE
    }

    name <- deparse(substitute(dds))
    contrast_name <- res_contrast_name(res)

    pheatmap(mat,
             annotation = annotation,
             clustering_distance_cols = "correlation",
             clustering_method = "ward.D2",
             main = paste(name, contrast_name, sep = " : "),
             scale = "row",
             show_rownames = show_rownames)
}
