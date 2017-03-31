#' Gene heatmap
#'
#' @author Michael Steinbaugh
#'
#' @import DESeq2
#' @import dplyr
#' @import pheatmap
#' @import SummarizedExperiment
#' @import tibble
#'
#' @param bcbio bcbio list object
#' @param dds \code{DESeqDataSet} object
#' @param de Plot differentially expressed (DE) genes. If \code{FALSE}, will
#'   plot a heatmap of all genes.
#'
#' @return Gene heatmap
#' @export
gene_heatmap <- function(bcbio,
                         dds,
                         alpha = 0.1,
                         lfc = 1,
                         de = TRUE) {
    check_bcbio_object(bcbio)
    if (class(dds)[1] != "DESeqDataSet") {
        stop("DESeqDataSet required")
    }

    name <- deparse(substitute(dds)) %>%
        gsub("_dds$", "", .)
    annotation <- import_metadata(bcbio) %>%
        .[, bcbio$intgroup]

    # Get a list of the DE genes
    res <- DESeq2::results(dds, alpha = alpha)

    # rlog transform the data
    rld <- DESeq2::rlog(dds)
    mat <- SummarizedExperiment::assay(rld)

    # Subset differentially expressed genes
    if (isTRUE(de)) {
        res <- res %>%
            as.data.frame %>%
            tibble::rownames_to_column("ensembl_gene_id") %>%
            dplyr::filter_(.dots = ~padj < alpha) %>%
            dplyr::filter_(.dots = ~log2FoldChange < -lfc |
                               log2FoldChange > lfc)
        # Subset the matrix
        mat <- mat[res$ensembl_gene_id, ]
    }

    pheatmap::pheatmap(mat,
                       annotation = annotation,
                       clustering_distance_cols = "correlation",
                       clustering_method = "ward.D2",
                       main = name,
                       scale = "row",
                       show_rownames = FALSE)
}
