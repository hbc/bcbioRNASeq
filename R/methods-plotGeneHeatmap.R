#' Gene Heatmap
#'
#' These functions facilitate heatmap plotting of a specified set of genes. By
#' default, row- and column-wise hierarchical clustering is performed using the
#' Ward method, but this behavior can be overrided by setting `cluster_rows` or
#' `cluster_cols` to `FALSE`. When column clustering is disabled, the columns
#' are sorted by the interesting groups (`interestingGroups`) specified in the
#' [bcbioRNADataSet] and then the sample names.
#'
#' @rdname plotGeneHeatmap
#' @name plotGeneHeatmap
#' @family Heatmaps
#' @author Michael Steinbaugh
#'
#' @inheritParams AllGenerics
#' @param annotationCol [data.frame] that specifies the annotations shown on the
#'   right side of the heatmap. Each row of this [data.frame] defines the
#'   features of the heatmap columns.
#' @param genes Character vector of specific gene identifiers to plot.
#' @param title *Optional*. Plot title.
#'
#' @seealso [pheatmap::pheatmap()].
#'
#' @return Graphical output only.
#'
#' @examples
#' data(bcb, dds, rld)
#'
#' # bcbioRNADataSet
#' plotGeneHeatmap(bcb)
#'
#' # Genes as Ensembl identifiers
#' genes <- counts(bcb)[1L:20L, ] %>% rownames
#' plotGeneHeatmap(bcb, genes = genes)
#'
#' # DESeqDataSet
#' plotGeneHeatmap(dds)
#'
#' # DESeqTransform
#' plotGeneHeatmap(rld)
NULL



# Constructors ====
.plotGeneHeatmap <- function(
    counts,
    genes = NULL,
    annotationCol = NULL,
    title = NULL,
    # Internal parameters
    scale = "row") {
    counts <- as.matrix(counts)

    # Check for identifier mismatch. Do this before zero count subsetting.
    if (!is.null(genes)) {
        if (!all(genes %in% rownames(counts))) {
            stop(paste(
                "Genes missing from counts matrix:",
                toString(setdiff(genes, rownames(counts)))),
                call. = FALSE)
        }
        counts <- counts %>%
            .[rownames(.) %in% genes, , drop = FALSE]
    }

    # Subset zero counts
    counts <- counts %>%
        .[rowSums(.) > 0L, , drop = FALSE]
    if (!is.matrix(counts) |
        nrow(counts) < 2L) {
        stop("Need at least 2 genes to cluster")
    }

    # Convert Ensembl gene identifiers to symbol names, if necessary
    if (nrow(counts) <= 100L) {
        showRownames <- TRUE
    } else {
        showRownames <- FALSE
    }
    if (isTRUE(showRownames)) {
        counts <- gene2symbol(counts)
    }

    if (!is.null(annotationCol)) {
        annotationCol <- annotationCol %>%
            as.data.frame %>%
            # Coerce annotation columns to factors
            rownames_to_column %>%
            mutate_all(factor) %>%
            column_to_rownames
        # Define colors for each annotation column
        annotationColors <- lapply(
            seq_along(dim(annotationCol)[[2L]]), function(a) {
                col <- annotationCol[[a]] %>%
                    levels
                colors <- annotationCol[[a]] %>%
                    levels %>%
                    length %>%
                    viridis
                names(colors) <- col
                colors
            }) %>%
            set_names(colnames(annotationCol))
    } else {
        annotationColors <- NULL
    }

    # pheatmap will error if `NULL` title is passed as `main`
    if (is.null(title)) {
        title <- ""
    }

    pheatmap(counts,
             annotation_col = annotationCol,
             annotation_colors = annotationColors,
             border_color = NA,
             color = inferno(256L),
             main = title,
             scale = scale,
             show_rownames = showRownames)
}



# Methods ====
#' @rdname plotGeneHeatmap
#' @export
setMethod("plotGeneHeatmap", "bcbioRNADataSet", function(
    object,
    genes = NULL,
    title = NULL) {
    counts <- counts(object, normalized = "rlog")
    annotationCol <- colData(object) %>%
        .[, metadata(object)[["interestingGroups"]], drop = FALSE]
    .plotGeneHeatmap(
        counts = counts,
        annotationCol = annotationCol,
        # User-defined
        genes = genes,
        title = title)
})



#' @rdname plotGeneHeatmap
#' @export
setMethod("plotGeneHeatmap", "DESeqDataSet", function(
    object,
    genes = NULL,
    annotationCol = NULL,
    title = NULL) {
    counts <- counts(object, normalized = TRUE)
    .plotGeneHeatmap(
        counts = counts,
        # User-defined
        genes = genes,
        annotationCol = annotationCol,
        title = title)
})



#' @rdname plotGeneHeatmap
#' @export
setMethod("plotGeneHeatmap", "DESeqTransform", function(
    object,
    genes = NULL,
    annotationCol = NULL,
    title = NULL) {
    counts <- assay(object)
    .plotGeneHeatmap(
        counts = counts,
        # User-defined
        genes = genes,
        annotationCol = annotationCol,
        title = title)
})



#' @rdname plotGeneHeatmap
#' @export
setMethod("plotGeneHeatmap", "matrix", function(
    object,
    genes = NULL,
    annotationCol = NULL,
    title = NULL) {
    .plotGeneHeatmap(
        counts = object,
        # User-defined
        genes = genes,
        annotationCol = annotationCol,
        title = title)
})
