#' Gene heatmap
#'
#' @author Michael Steinbaugh
#'
#' @param counts Counts matrix
#' @param deg Optional data frame of differential expressed genes (DEG) produced
#'   by DESeq2. If \code{NULL}, will output a heatmap of all genes.
#' @param metadata Metadata data frame
#' @param factor Factor list
gene_heatmap <- function(counts,
                         deg = NULL,
                         metadata,
                         factor = "intgroup") {
    metadata <- as.data.frame(metadata)

    if (!is.data.frame(metadata)) {
        stop("A metadata data frame is required.")
    }
    if (!is.character(factor)) {
        stop("A factor group character vector is required.")
    }
    if (class(counts)[1] == "DESeqTransform") {
        counts <- SummarizedExperiment::assay(counts)
    }
    if (!is.matrix(counts)) {
        stop("A counts matrix is required.")
    }

    if (!is.null(deg)) {
        if (!is.data.frame(deg)) {
            stop("DEG must be input as a data frame.")
        }
        name <- deparse(substitute(deg))
        counts <- counts[deg$ensembl_gene_id, ]
    } else {
        name <- "all genes"
        # Error in hclust(d, method = method) :
        #     NA/NaN/Inf in foreign function call (arg 11)
    }

    pheatmap::pheatmap(counts,
                       annotation = metadata[, factor],
                       clustering_distance_cols = "correlation",
                       clustering_method = "ward.D2",
                       main = name,
                       scale = "row",
                       show_rownames = FALSE)
}
