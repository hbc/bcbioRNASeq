#' Differential expression heatmap
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
#' @param contrast Contrast
#' @param alpha Alpha level cutoff
#' @param lfc Log fold change ratio (base 2) cutoff
#'
#' @return Gene heatmap
#' @export
de_heatmap <- function(bcbio,
                       dds,
                       contrast,
                       alpha = 0.1,
                       lfc = 0) {
    check_bcbio_object(bcbio)
    if (class(dds)[1] != "DESeqDataSet") {
        stop("DESeqDataSet required")
    }

    name <- deparse(substitute(dds)) %>%
        gsub("_dds$", "", .)
    annotation <- import_metadata(bcbio) %>%
        .[, bcbio$intgroup]

    # Obtain the DE genes
    res <- DESeq2::results(dds,
                           contrast = contrast,
                           alpha = alpha) %>%
        as.data.frame %>%
        tibble::rownames_to_column("ensembl_gene_id") %>%
        dplyr::filter_(.dots = ~padj < alpha) %>%
        dplyr::filter_(.dots = ~log2FoldChange < -lfc |
                           log2FoldChange > lfc)

    # rlog transform the data
    rld <- DESeq2::rlog(dds)
    mat <- rld %>%
        SummarizedExperiment::assay(.) %>%
        .[res$ensembl_gene_id, ]

    pheatmap::pheatmap(mat,
                       annotation = annotation,
                       clustering_distance_cols = "correlation",
                       clustering_method = "ward.D2",
                       main = name,
                       scale = "row",
                       show_rownames = FALSE)
}
