#' Plot Heatmap
#'
#' These functions facilitate heatmap plotting of a specified set of genes. By
#' default, row- and column-wise hierarchical clustering is performed using the
#' Ward method, but this behavior can be overrided by setting `cluster_rows` or
#' `cluster_cols` to `FALSE`. When column clustering is disabled, the columns
#' are sorted by the interesting groups (`interestingGroups`) specified in the
#' [bcbioRNASeq] and then the sample names.
#'
#' @rdname plotHeatmap
#' @name plotHeatmap
#' @family Heatmaps
#' @author Michael Steinbaugh
#'
#' @inheritParams AllGenerics
#' @inheritParams basejump::gene2symbol
#' @inheritParams counts
#'
#' @param genes *Optional*. Character vector of gene identifiers to plot. These
#'   must be the stable identifiers (e.g. ENSG00000000003) used on Ensembl and
#'   not the gene symbols.
#' @param gene2symbol Apply gene identifier to symbol mappings. If set `TRUE`,
#'   the function will attempt to automatically map gene identifiers to symbols
#'   from Ensembl using [basejump::annotable()]. If set `FALSE`/`NULL`, then
#'   gene2symbol mapping will be disabled. This is useful when working with a
#'   poorly annotated genome. Alternatively, a gene2symbol [data.frame] can be
#'   passed in, and must contain the columns `ensgene` and `symbol`. then the
#'   Ensembl gene identifiers will be labeled in place of gene symbols.
#' @param annotationCol *Optional*. [data.frame] that specifies the annotations
#'   shown on the right side of the heatmap. Each row of this [data.frame]
#'   defines the features of the heatmap columns.
#' @param scale Character indicating if the values should be centered and scaled
#'   in either the row direction or the column direction, or none. Corresponding
#'   values are "row", "column" and "none".
#' @param color Colors to use for plot. Defaults to the [viridis::viridis()]
#'   palette.
#' @param legendColor Colors to use for legend labels. Defaults to the
#'   [viridis::viridis()] palette.
#' @param title *Optional*. Plot title.
#' @param ... Passthrough arguments to [pheatmap::pheatmap()].
#'
#' @seealso [pheatmap::pheatmap()].
#'
#' @return Graphical output only.
#'
#' @examples
#' bcb <- examples[["bcb"]]
#'
#' # Use Ensembl identifiers to define genes
#' genes <- counts(bcb)[1:20, ] %>% rownames()
#' plotHeatmap(bcb, genes = genes)
#'
#' # Use inferno color palette
#' plotHeatmap(
#'     bcb,
#'     genes = genes,
#'     color = viridis::inferno(256),
#'     legendColor = viridis::inferno)
#'
#' # Transcriptome heatmap with default pheatmap colors
#' plotHeatmap(
#'     bcb,
#'     color = NULL,
#'     legendColor = NULL)
#'
#' # DESeqDataSet
#' dds <- examples[["dds"]]
#' plotHeatmap(dds)
#'
#' # DESeqTransform
#' rld <- examples[["rld"]]
#' plotHeatmap(rld)
NULL



# Constructors ====
#' @importFrom basejump checkGene2symbol gene2symbol
#' @importFrom dplyr mutate_all
#' @importFrom pheatmap pheatmap
#' @importFrom stats setNames
#' @importFrom tibble column_to_rownames rownames_to_column
#' @importFrom viridis viridis
.plotHeatmap <- function(
    object,
    genes = NULL,
    gene2symbol = TRUE,
    annotationCol = NULL,
    scale = "row",
    color = viridis::viridis(256),
    legendColor = viridis::viridis,
    title = NULL,
    quiet = FALSE,
    ...) {
    counts <- object
    if (!is.matrix(counts)) return(NULL)

    # Check for missing genes
    if (!is.null(genes)) {
        if (!all(genes %in% rownames(counts))) {
            stop(paste(
                "Genes missing from counts matrix:",
                toString(setdiff(genes, rownames(counts)))),
                call. = FALSE)
        }
        counts <- counts %>%
            .[rownames(.) %in% genes, , drop = FALSE]
    } else {
        # Remove zero counts
        counts <- counts %>%
            .[rowSums(.) > 0, , drop = FALSE]
    }

    if (nrow(counts) < 2) {
        stop("Need at least 2 genes to plot heatmap", call. = FALSE)
    }

    # Convert Ensembl gene identifiers to symbol names, if necessary
    if (nrow(counts) <= 100) {
        showRownames <- TRUE
    } else {
        showRownames <- FALSE
    }
    if (isTRUE(showRownames)) {
        if (isTRUE(gene2symbol)) {
            counts <- gene2symbol(counts, quiet = quiet)
        } else if (is.data.frame(gene2symbol)) {
            checkGene2symbol(gene2symbol)
            # Remap the rownames to use the gene symbols
            remap <- gene2symbol %>%
                .[rownames(counts), "symbol"] %>%
                make.names(unique = TRUE)
            rownames(counts) <- remap
        }
    }

    # Prepare the annotation columns
    if (!is.null(annotationCol)) {
        annotationCol <- annotationCol %>%
            as.data.frame() %>%
            rownames_to_column() %>%
            mutate_all(factor) %>%
            column_to_rownames()
    }

    # Define colors for each annotation column, if desired
    if (!is.null(annotationCol) & !is.null(legendColor)) {
        annotationColors <- lapply(
            seq_along(colnames(annotationCol)), function(a) {
                col <- annotationCol[[a]] %>%
                    levels()
                colors <- annotationCol[[a]] %>%
                    levels() %>%
                    length() %>%
                    legendColor
                names(colors) <- col
                colors
            }) %>%
            setNames(colnames(annotationCol))
    } else {
        annotationColors <- NULL
    }

    # pheatmap will error if `NULL` title is passed as `main`
    if (is.null(title)) {
        title <- ""
    }

    # If `color = NULL`, use the pheatmap default
    if (is.null(color)) {
        color <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
    }

    pheatmap(
        mat = counts,
        annotation_col = annotationCol,
        annotation_colors = annotationColors,
        border_color = NA,
        color = color,
        main = title,
        scale = scale,
        show_rownames = showRownames,
        ...)
}



# Methods ====
#' @rdname plotHeatmap
#' @importFrom viridis viridis
#' @export
setMethod(
    "plotHeatmap",
    signature("bcbioRNASeq"),
    function(
        object,
        normalized = "rlog",
        genes = NULL,
        scale = "row",
        color = viridis::viridis(256),
        legendColor = viridis::viridis,
        title = NULL,
        quiet = FALSE,
        ...) {
        counts <- counts(object, normalized = normalized)
        gene2symbol <- gene2symbol(object)
        interestingGroups <- interestingGroups(object)
        annotationCol <- colData(object) %>%
            .[, interestingGroups, drop = FALSE]
        .plotHeatmap(
            object = counts,
            genes = genes,
            gene2symbol = gene2symbol,
            annotationCol = annotationCol,
            scale = scale,
            color = color,
            legendColor = legendColor,
            title = title,
            quiet = quiet,
            ...)
    })



#' @rdname plotHeatmap
#' @importFrom viridis viridis
#' @export
setMethod(
    "plotHeatmap",
    signature("DESeqDataSet"),
    function(
        object,
        normalized = TRUE,
        genes = NULL,
        gene2symbol = TRUE,
        annotationCol = NULL,
        scale = "row",
        color = viridis::viridis(256),
        legendColor = viridis::viridis,
        title = NULL,
        quiet = FALSE,
        ...) {
        counts <- counts(object, normalized = normalized)
        .plotHeatmap(
            object = counts,
            genes = genes,
            gene2symbol = gene2symbol,
            annotationCol = annotationCol,
            scale = scale,
            color = color,
            legendColor = legendColor,
            title = title,
            quiet = quiet,
            ...)
    })



#' @rdname plotHeatmap
#' @importFrom viridis viridis
#' @export
setMethod(
    "plotHeatmap",
    signature("DESeqTransform"),
    function(
        object,
        genes = NULL,
        gene2symbol = TRUE,
        annotationCol = NULL,
        scale = "row",
        color = viridis::viridis(256),
        legendColor = viridis::viridis,
        title = NULL,
        quiet = FALSE,
        ...) {
        counts <- assay(object)
        .plotHeatmap(
            object = counts,
            genes = genes,
            gene2symbol = gene2symbol,
            annotationCol = annotationCol,
            scale = scale,
            color = color,
            legendColor = legendColor,
            title = title,
            quiet = quiet,
            ...)
    })



#' @rdname plotHeatmap
#' @importFrom viridis viridis
#' @export
setMethod(
    "plotHeatmap",
    signature("matrix"),
    .plotHeatmap)
