#' Plot Individual Genes
#'
#' @rdname plotGene
#' @name plotGene
#' @family Quality Control Plots
#' @author Michael Steinbaugh
#'
#' @importFrom bcbioBase plotGene
#'
#' @inherit plotTotalReads
#'
#' @param genes *Optional*. Gene identifiers (rownames) to plot. These must be
#'   the stable identifiers (e.g. ENSG00000000003) used on Ensembl and not the
#'   gene symbols.
#' @param gene2symbol Apply gene identifier to symbol mappings. A gene2symbol
#'   [data.frame] can be passed in, and must contain the columns `ensgene` and
#'   `symbol`. then the Ensembl gene identifiers will be labeled in place of
#'   gene symbols.
#' @param metadata Sample metadata [data.frame].
#' @param stackReplicates Stack replicate points into a single tick on the
#'   sample axis.
#' @param color Desired ggplot color scale. Defaults to
#'   [viridis::scale_color_viridis()]. Must supply discrete values. When set to
#'   `NULL`, the default ggplot2 color palette will be used. If manual color
#'   definitions are desired, we recommend using
#'   [ggplot2::scale_color_manual()].
#' @param countsAxisLabel Text label of counts axis.
#' @param return Desired return type: `grid`, `list`, `markdown`.
#' @param headerLevel R Markdown header level. Only applies when
#'   `return = "markdown"`.
#'
#' @return
#' - `grid`: [cowplot::plot_grid()] graphical output.
#' - `list`: Plot list, containing per gene [ggplot] objects.
#' - `markdown`: Tabset R Markdown output, tabbed per gene.
#'
#' @seealso [DESeq2::plotCounts()].
#'
#' @examples
#' load(system.file(
#'     file.path("extdata", "bcb.rda"),
#'     package = "bcbioRNASeq"))
#'
#' # Gene identifiers
#' genes <- rownames(bcb)[1:4]
#' print(genes)
#' plotGene(bcb, genes = genes)
#'
#' # Default ggplot2 color palette
#' plotGene(
#'     bcb,
#'     genes = genes,
#'     interestingGroups = "sampleName",
#'     color = NULL)
NULL



# Constructors =================================================================
#' Plot Gene Constructor
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @importFrom bcbioBase annotable detectOrganism uniteInterestingGroups
#' @importFrom cowplot plot_grid
#' @importFrom ggplot2 aes_string element_text expand_limits geom_point ggplot
#'   guides labs theme
#' @importFrom pbapply pblapply
#' @importFrom tibble tibble
#' @importFrom viridis scale_color_viridis
#'
#' @param genes Gene identifiers, as a named character vector. `ensgene` is
#'   provided as the value and `symbol` as the name. This gets defined in the S4
#'   method (see below).
#'
#' @return [ggplot].
.plotGene <- function(
    object,
    genes,
    metadata,
    gene2symbol = NULL,
    interestingGroups = "sampleName",
    stackReplicates = TRUE,
    color = viridis::scale_color_viridis(discrete = TRUE),
    countsAxisLabel = "counts",
    return = "grid",
    headerLevel = 2L) {
    if (!(is.data.frame(gene2symbol) || is.null(gene2symbol))) {
        abort("`gene2symbol` must be a data.frame or NULL")
    }

    # Gene to symbol mappings
    if (is.data.frame(gene2symbol)) {
        checkGene2symbol(gene2symbol)
        if (!all(genes %in% gene2symbol[["ensgene"]])) {
            abort(paste(
                "Genes missing in gene2symbol:",
                setdiff(genes, gene2symbol[["ensgene"]])
            ))
        }
        match <- match(x = genes, table = gene2symbol[["ensgene"]])
        gene2symbol <- gene2symbol[match, , drop = FALSE]
        genes <- gene2symbol[["ensgene"]]
        names(genes) <- gene2symbol[["symbol"]]
    }

    # Prepare interesting groups column
    metadata <- metadata %>%
        as.data.frame() %>%
        uniteInterestingGroups(interestingGroups)

    plots <- lapply(seq_along(genes), function(a) {
        ensgene <- genes[[a]]
        symbol <- names(genes)[[a]]
        if (!is.null(symbol)) {
            title <- symbol
            subtitle <- ensgene
        } else {
            title <- ensgene
            subtitle <- NULL
        }
        if (isTRUE(stackReplicates)) {
            x <- metadata[["interestingGroups"]]
            xlab <- paste(interestingGroups, collapse = ":")
        } else {
            x <- colnames(object)
            xlab <- "sample"
        }
        data <- tibble(
            x = x,
            y = object[ensgene, , drop = TRUE],
            interestingGroups = metadata[["interestingGroups"]])
        p <- ggplot(
            data,
            mapping = aes_string(
                x = "x",
                y = "y",
                color = "interestingGroups")
        ) +
            geom_point(size = 4L) +
            theme(axis.text.x = element_text(angle = 90L)) +
            labs(
                title = title,
                subtitle = subtitle,
                x = xlab,
                y = countsAxisLabel,
                color = paste(interestingGroups, collapse = ":\n")
            ) +
            expand_limits(y = 0L)
        if (is(color, "ScaleDiscrete")) {
            p <- p + color
        }
        if (identical(interestingGroups, "sampleName")) {
            p <- p + guides(color = FALSE)
        }
        p
    })
    names(plots) <- genes

    # Return
    validReturn <- c("grid", "list", "markdown")
    if (return == "grid") {
        plot_grid(plotlist = plots, labels = NULL)
    } else if (return == "list") {
        plots
    } else if (return == "markdown") {
        pblapply(seq_along(plots), function(a) {
            if (is.numeric(headerLevel)) {
                ensgene <- genes[[a]]
                symbol <- names(genes)[[a]]
                if (is.null(symbol)) {
                    symbol <- ensgene
                }
                mdHeader(symbol, level = headerLevel, asis = TRUE)
            }
            show(plots[[a]])
        }) %>%
            invisible()
    } else {
        abort(paste("`return` must contain:", toString(validReturn)))
    }
}



# Methods ======================================================================
#' @rdname plotGene
#' @importFrom viridis scale_color_viridis
#' @param normalized Normalization method. Supports `tpm` (**default**), `tmm`,
#'   `rlog`, or `vst`.
#' @export
setMethod(
    "plotGene",
    signature("bcbioRNASeq"),
    function(
        object,
        genes,
        normalized = "tpm",
        log2 = TRUE,
        interestingGroups,
        stackReplicates = TRUE,
        color = viridis::scale_color_viridis(discrete = TRUE),
        return = "grid",
        headerLevel = 2L) {
        if (identical(normalized, FALSE)) {
            abort(paste(
                "Raw counts are not library size adjusted and therefore",
                "not recommended to use for gene plots"
            ))
        }
        # Passthrough: stackReplicates, color, return, headerLevel
        if (missing(interestingGroups)) {
            interestingGroups <- bcbioBase::interestingGroups(object)
        }
        counts <- counts(object, normalized = normalized)
        if (isTRUE(normalized)) {
            countsAxisLabel <- "normalized counts"
            if (isTRUE(log2)) {
                counts <- log2(counts + 1L)
                countsAxisLabel <- c("log2", countsAxisLabel)
            }
        } else {
            countsAxisLabel <- normalized
        }
        gene2symbol <- gene2symbol(object)
        metadata <- sampleMetadata(object)
        .plotGene(
            object = counts,
            genes = genes,
            gene2symbol = gene2symbol,
            metadata = metadata,
            interestingGroups = interestingGroups,
            color = color,
            countsAxisLabel = countsAxisLabel,
            return = return,
            headerLevel = headerLevel)
    })



#' @rdname plotGene
#' @param log2 Apply `log2(x + 1)` transformation to normalized counts.
#' @export
setMethod(
    "plotGene",
    signature("DESeqDataSet"),
    function(
        object,
        genes,
        log2 = TRUE,
        gene2symbol = NULL,
        interestingGroups = "sampleName",
        stackReplicates = TRUE,
        color = viridis::scale_color_viridis(discrete = TRUE),
        return = "grid",
        headerLevel = 2L) {
        # Passthrough: genes, gene2symbol, interestingGroups, stackReplicates,
        # color, return, headerLevel
        counts <- counts(object, normalized = TRUE)
        countsAxisLabel <- "normalized counts"
        if (isTRUE(log2)) {
            counts <- log2(counts + 1L)
            countsAxisLabel <- paste("log2", countsAxisLabel)
        }
        metadata <- colData(object)
        .plotGene(
            object = counts,
            genes = genes,
            gene2symbol = gene2symbol,
            metadata = metadata,
            interestingGroups = interestingGroups,
            stackReplicates = stackReplicates,
            color = color,
            countsAxisLabel = countsAxisLabel,
            return = return,
            headerLevel = headerLevel)
    })



#' @rdname plotGene
#' @export
setMethod(
    "plotGene",
    signature("DESeqTransform"),
    function(
        object,
        genes,
        gene2symbol = NULL,
        interestingGroups = "sampleName",
        stackReplicates = TRUE,
        color = viridis::scale_color_viridis(discrete = TRUE),
        return = "grid",
        headerLevel = 2L) {
        counts <- assay(object)
        metadata <- colData(object)
        # Passthrough: genes, gene2symbol, interestingGroups, stackReplicates,
        # color, return, headerLevel
        if ("rlogIntercept" %in% colnames(mcols(object))) {
            countsAxisLabel <- "rlog"
        } else {
            countsAxisLabel <- "vst"
        }
        .plotGene(
            object = counts,
            genes = genes,
            gene2symbol = gene2symbol,
            metadata = metadata,
            interestingGroups = interestingGroups,
            stackReplicates = stackReplicates,
            color = color,
            countsAxisLabel = countsAxisLabel,
            return = return,
            headerLevel = headerLevel)
    })



#' @rdname plotGene
#' @export
setMethod(
    "plotGene",
    signature("matrix"),
    .plotGene)
