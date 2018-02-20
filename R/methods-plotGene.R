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
#' @param genes Gene identifiers (rownames) to plot. These must be the stable
#'   identifiers (e.g. ENSG00000000003) used on Ensembl and not the gene
#'   symbols.
#' @param gene2symbol Apply gene identifier to symbol mappings. A gene2symbol
#'   [data.frame] can be passed in, and must contain the columns `ensgene` and
#'   `symbol`. then the Ensembl gene identifiers will be labeled in place of
#'   gene symbols.
#' @param metadata Sample metadata [data.frame].
#' @param normalized Normalization method. Supports `tpm` (**default**), `tmm`,
#'   `rlog`, or `vst`.
#' @param stackReplicates Stack replicate points into a single tick on the
#'   sample axis.
#' @param color Desired ggplot color scale. Defaults to
#'   [scale_color_viridis()]. Must supply discrete values. When set to
#'   `NULL`, the default ggplot2 color palette will be used. If manual color
#'   definitions are desired, we recommend using [scale_color_manual()].
#' @param countsAxisLabel Text label of counts axis.
#' @param return Desired return type: `grid`, `wide`, `list`, or `markdown`.
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
#' load(system.file("extdata/bcb.rda", package = "bcbioRNASeq"))
#' load(system.file("extdata/dds.rda", package = "bcbioRNASeq"))
#' load(system.file("extdata/rld.rda", package = "bcbioRNASeq"))
#'
#' # Gene identifiers
#' genes <- rownames(bcb)[1:4]
#'
#' # bcbioRNASeq
#' plotGene(bcb, genes = genes)
#' plotGene(
#'     bcb,
#'     genes = genes,
#'     interestingGroups = "sampleName",
#'     color = NULL)
#'
#' # DESeqDataSet
#' plotGene(dds, genes = genes, interestingGroups = "group")
#'
#' # DESeqTransform
#' plotGene(rld, genes = genes, interestingGroups = "group")
NULL



# Constructors =================================================================
#' @importFrom basejump annotable detectOrganism markdownHeader
#' @importFrom bcbioBase uniteInterestingGroups
#' @importFrom BiocParallel bplapply
#' @importFrom cowplot plot_grid
#' @importFrom ggplot2 aes_string element_text expand_limits geom_point ggplot
#'   guides labs theme
#' @importFrom tibble tibble
.plotGene <- function(
    object,
    genes,
    gene2symbol = NULL,
    metadata,
    interestingGroups = "sampleName",
    stackReplicates = TRUE,
    color = scale_color_viridis(discrete = TRUE),
    countsAxisLabel = "counts",
    return = "grid",
    headerLevel = 2L) {
    assert_is_character(genes)
    assert_is_data.frame(metadata)
    assertFormalGene2symbol(object, genes, gene2symbol)
    assert_is_data.frame(metadata)
    assertFormalInterestingGroups(metadata, interestingGroups)
    assert_is_a_bool(stackReplicates)
    assertIsColorScaleDiscreteOrNULL(color)
    assert_is_a_string(countsAxisLabel)
    assert_is_a_string(return)
    assert_is_subset(return, c("grid", "list", "markdown", "wide"))
    assertIsAHeaderLevel(headerLevel)

    # Gene to symbol mappings
    if (is.data.frame(gene2symbol)) {
        match <- match(x = genes, table = gene2symbol[["ensgene"]])
        gene2symbol <- gene2symbol[match, , drop = FALSE]
        genes <- gene2symbol[["ensgene"]]
        names(genes) <- gene2symbol[["symbol"]]
    }

    # Prepare interesting groups column
    metadata <- uniteInterestingGroups(metadata, interestingGroups)

    plots <- bplapply(seq_along(genes), function(a) {
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
            data = data,
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
    if (return == "grid") {
        plot_grid(plotlist = plots, labels = NULL)
    } else if (return == "list") {
        plots
    } else if (return == "markdown") {
        bplapply(seq_along(plots), function(a) {
            if (is.numeric(headerLevel)) {
                ensgene <- genes[[a]]
                symbol <- names(genes)[[a]]
                if (is.null(symbol)) {
                    symbol <- ensgene
                }
                markdownHeader(symbol, level = headerLevel, asis = TRUE)
            }
            show(plots[[a]])
        }) %>%
            invisible()
    }
}



.plotGene.bcbioRNASeq <- function(  # nolint
    object,
    genes,
    normalized = "rlog",
    interestingGroups,
    stackReplicates = TRUE,
    color = scale_color_viridis(discrete = TRUE),
    return = "grid",
    headerLevel = 2L) {
    # Passthrough: stackReplicates, color, return, headerLevel
    assert_is_a_string(normalized)
    if (missing(interestingGroups)) {
        interestingGroups <- bcbioBase::interestingGroups(object)
    }
    .plotGene(
        object = counts(object, normalized = normalized),
        genes = genes,
        gene2symbol = gene2symbol(object),
        metadata = sampleMetadata(object),
        interestingGroups = interestingGroups,
        color = color,
        countsAxisLabel = normalized,
        return = return,
        headerLevel = headerLevel)
}



.plotGene.DESeqDataSet <- function(  # nolint
    object,
    genes,
    gene2symbol = NULL,
    interestingGroups = "sampleName",
    stackReplicates = TRUE,
    color = scale_color_viridis(discrete = TRUE),
    return = "grid",
    headerLevel = 2L) {
    # Passthrough: genes, gene2symbol, interestingGroups, stackReplicates,
    # color, return, headerLevel
    .plotGene(
        object = log2(counts(object, normalized = TRUE) + 1L),
        genes = genes,
        gene2symbol = gene2symbol,
        metadata = sampleMetadata(object),
        interestingGroups = interestingGroups,
        stackReplicates = stackReplicates,
        color = color,
        countsAxisLabel = "log2 normalized counts",
        return = return,
        headerLevel = headerLevel)
}



.plotGene.DESeqTransform <- function(  # nolint
    object,
    genes,
    gene2symbol = NULL,
    interestingGroups = "sampleName",
    stackReplicates = TRUE,
    color = scale_color_viridis(discrete = TRUE),
    return = "grid",
    headerLevel = 2L) {
    # Passthrough: genes, gene2symbol, interestingGroups, stackReplicates,
    # color, return, headerLevel
    if ("rlogIntercept" %in% colnames(mcols(object))) {
        countsAxisLabel <- "rlog"
    } else {
        countsAxisLabel <- "vst"
    }
    .plotGene(
        object = assay(object),
        genes = genes,
        gene2symbol = gene2symbol,
        metadata = sampleMetadata(object),
        interestingGroups = interestingGroups,
        stackReplicates = stackReplicates,
        color = color,
        countsAxisLabel = countsAxisLabel,
        return = return,
        headerLevel = headerLevel)
}



# Methods ======================================================================
#' @rdname plotGene
#' @export
setMethod(
    "plotGene",
    signature("bcbioRNASeq"),
    .plotGene.bcbioRNASeq)



#' @rdname plotGene
#' @export
setMethod(
    "plotGene",
    signature("DESeqDataSet"),
    .plotGene.DESeqDataSet)



#' @rdname plotGene
#' @export
setMethod(
    "plotGene",
    signature("DESeqTransform"),
    .plotGene.DESeqTransform)



#' @rdname plotGene
#' @export
setMethod(
    "plotGene",
    signature("matrix"),
    .plotGene)
