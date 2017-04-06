#' Differentially expressed gene (DEG) heatmap
#'
#' @author Michael Steinbaugh
#'
#' @import DESeq2
#' @import dplyr
#' @import pheatmap
#' @import SummarizedExperiment
#' @import tibble
#' @importFrom S4Vectors mcols
#'
#' @param bcbio bcbio list object
#' @param dds \code{DESeqDataSet} object
#' @param contrast Contrast
#' @param alpha Alpha level cutoff
#' @param lfc Log fold change ratio (base 2) cutoff
#'
#' @return Gene heatmap
#' @export
deg_heatmap <- function(
    bcbio,
    dds,
    contrast,
    alpha = 0.1,
    lfc = 0) {
    check_bcbio(bcbio)
    check_dds(dds)

    # rlog transform the counts
    rld <- DESeq2::rlog(dds)

    # Annotation metadata
    annotation <- dds %>%
        SummarizedExperiment::colData(.) %>%
        as.data.frame %>%
        .[, bcbio$intgroup]

    # DE genes, used for subsetting the rlog counts matrix
    res <- DESeq2::results(dds,
                           contrast = contrast,
                           alpha = alpha)

    # Output results to dataframe and subset by alpha and lfc cutoffs
    res_df <- res %>%
        as.data.frame %>%
        tibble::rownames_to_column("ensembl_gene_id") %>%
        dplyr::filter_(.dots = ~padj < alpha) %>%
        dplyr::filter_(.dots = ~log2FoldChange < -lfc | log2FoldChange > lfc)

    # Subset the transformed count matrix
    mat <- rld %>%
        SummarizedExperiment::assay(.) %>%
        .[res_df$ensembl_gene_id, ]

    if (nrow(mat) <= 100) {
        show_rownames <- TRUE
        # Change rownames to readable external gene names
        gene_names <- ensembl_annotations(bcbio, values = rownames(mat))
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

    pheatmap::pheatmap(mat,
                       annotation = annotation,
                       clustering_distance_cols = "correlation",
                       clustering_method = "ward.D2",
                       main = paste(name, contrast_name, sep = " : "),
                       scale = "row",
                       show_rownames = show_rownames)
}
