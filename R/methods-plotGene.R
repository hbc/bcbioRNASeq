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
#' @param colData Sample metadata [data.frame].
#' @param normalized Normalization method. Supports `tpm` (**default**), `tmm`,
#'   `rlog`, or `vst`.
#' @param countsAxisLabel Text label of counts axis.
#' @param medianLine Include median line for each group. Disabled when samples
#'   are grouped by `sampleName`.
#' @param color Desired ggplot color scale. Defaults to
#'   [scale_color_viridis()]. Must supply discrete values. When set to
#'   `NULL`, the default ggplot2 color palette will be used. If manual color
#'   definitions are desired, we recommend using [scale_color_manual()].
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
#' genes <- head(rownames(bcb), 4L)
#'
#' # bcbioRNASeq ====
#' plotGene(bcb, genes = genes, return = "grid")
#' plotGene(bcb, genes = genes, return = "wide")
#'
#' # DESeqDataSet ====
#' plotGene(dds, genes = genes, interestingGroups = "group")
#'
#' # DESeqTransform ====
#' plotGene(rld, genes = genes, interestingGroups = "group")
NULL



# Constructors =================================================================
#' @importFrom basejump markdownPlotlist
#' @importFrom bcbioBase uniteInterestingGroups
#' @importFrom cowplot plot_grid
.plotGene <- function(
    object,
    genes,
    gene2symbol = NULL,
    colData,
    interestingGroups = "sampleName",
    countsAxisLabel = "counts",
    medianLine = TRUE,
    color = scale_color_viridis(discrete = TRUE),
    return = "grid",
    headerLevel = 2L) {
    assert_is_matrix(object)
    assert_has_dimnames(object)
    assert_is_character(genes)
    assert_is_subset(genes, rownames(object))
    object <- object[genes, , drop = FALSE]
    assertFormalGene2symbol(object, genes, gene2symbol)
    assertFormalAnnotationCol(object, colData)
    assertFormalInterestingGroups(colData, interestingGroups)
    assert_is_a_string(countsAxisLabel)
    assert_is_a_bool(medianLine)
    assertIsColorScaleDiscreteOrNULL(color)
    assert_is_a_string(return)
    assert_is_subset(return, c("grid", "list", "markdown", "wide"))
    assertIsAHeaderLevel(headerLevel)

    # Add `interestingGroups` column to colData
    colData <- uniteInterestingGroups(colData, interestingGroups)

    if (return != "wide") {
        plotlist <- .plotGeneList(
            object = object,
            genes = genes,
            gene2symbol = gene2symbol,
            colData = colData,
            interestingGroups = interestingGroups,
            countsAxisLabel = countsAxisLabel,
            medianLine = medianLine,
            color = color)
    }

    if (return == "grid") {
        if (length(plotlist) > 1L) {
            labels <- "AUTO"
        } else {
            labels <- NULL
        }
        plot_grid(plotlist = plotlist, labels = labels)
    } else if (return == "wide") {
        .plotGeneWide(
            object = object,
            genes = genes,
            gene2symbol = gene2symbol,
            colData = colData,
            interestingGroups = interestingGroups,
            countsAxisLabel = countsAxisLabel,
            medianLine = medianLine,
            color = color)
    } else if (return == "list") {
        plotlist
    } else if (return == "markdown") {
        markdownPlotlist(plotlist, headerLevel = headerLevel)
    }
}



#' @importFrom BiocParallel bpmapply
#' @importFrom ggplot2 aes_string element_text expand_limits geom_point ggplot
#'   guides labs theme
#' @importFrom tibble tibble
.plotGeneList <- function(
    object,
    genes,
    gene2symbol = NULL,
    colData,
    interestingGroups = "sampleName",
    countsAxisLabel = "counts",
    medianLine = TRUE,
    color = scale_color_viridis(discrete = TRUE)) {
    assert_is_subset(genes, rownames(object))
    object <- object[genes, , drop = FALSE]

    # Gene to symbol mappings
    if (is.data.frame(gene2symbol)) {
        assert_is_subset(genes, gene2symbol[["ensgene"]])
        match <- match(x = genes, table = gene2symbol[["ensgene"]])
        symbols <- gene2symbol[match, "symbol", drop = TRUE]
        rownames(object) <- symbols
        genes <- symbols
    }

    bpmapply(
        gene = genes,
        MoreArgs = list(
            object = object,
            colData = colData,
            interestingGroups = interestingGroups,
            countsAxisLabel = countsAxisLabel,
            medianLine = medianLine,
            color = color
        ),
        FUN = function(
            object,
            gene,
            colData,
            interestingGroups,
            countsAxisLabel,
            medianLine,
            color
        ) {
            data <- tibble(
                x = colData[["interestingGroups"]],
                y = object[gene, , drop = TRUE],
                interestingGroups = colData[["interestingGroups"]])

            p <- ggplot(
                data = data,
                mapping = aes_string(
                    x = "x",
                    y = "y",
                    color = "interestingGroups")
            ) +
                genePoint() +
                theme(axis.text.x = element_text(angle = 90L)) +
                labs(
                    title = gene,
                    x = NULL,
                    y = countsAxisLabel,
                    color = paste(interestingGroups, collapse = ":\n")
                ) +
                expand_limits(y = 0L) +
                theme(legend.position = "none")

            if (
                isTRUE(medianLine) &&
                !identical(interestingGroups, "sampleName")
            ) {
                p <- p + geneMedianLine
            }

            if (is(color, "ScaleDiscrete")) {
                p <- p + color
            }

            p
        },
        SIMPLIFY = FALSE,
        USE.NAMES = TRUE)
}



#' @importFrom dplyr arrange group_by left_join
#' @importFrom ggplot2 aes_string element_text expand_limits geom_point ggplot
#'   guides labs theme
#' @importFrom magrittr set_colnames
#' @importFrom rlang !! !!! sym syms
#' @importFrom reshape2 melt
#' @importFrom tibble as_tibble rownames_to_column
.plotGeneWide <- function(
    object,
    genes,
    gene2symbol = NULL,
    colData,
    interestingGroups = "sampleName",
    countsAxisLabel = "counts",
    medianLine = TRUE,
    color = scale_color_viridis(discrete = TRUE),
    title = NULL) {
    # Gene to symbol mappings
    if (is.data.frame(gene2symbol)) {
        assert_is_subset(genes, gene2symbol[["ensgene"]])
        match <- match(x = genes, table = gene2symbol[["ensgene"]])
        symbols <- gene2symbol[match, "symbol", drop = TRUE]
        rownames(object) <- symbols
    }

    # Melt counts into long format
    data <- object %>%
        as.data.frame() %>%
        rownames_to_column() %>%
        melt(id = 1L) %>%
        as_tibble() %>%
        set_colnames(c("gene", "sampleID", "counts")) %>%
        arrange(!!!syms(c("gene", "sampleID"))) %>%
        group_by(!!sym("gene")) %>%
        left_join(colData, by = "sampleID")

    p <- ggplot(
        data = data,
        mapping = aes_string(
            x = "gene",
            y = "counts",
            color = "interestingGroups")
    ) +
        genePoint() +
        theme(axis.text.x = element_text(angle = 90L, hjust = 1L)) +
        labs(
            title = title,
            x = NULL,
            y = countsAxisLabel,
            color = paste(interestingGroups, collapse = ":\n")
        ) +
        expand_limits(y = 0L)

    if (isTRUE(medianLine) && !identical(interestingGroups, "sampleName")) {
        p <- p + geneMedianLine
    }

    if (is(color, "ScaleDiscrete")) {
        p <- p + color
    }

    if (identical(interestingGroups, "sampleName")) {
        p <- p + guides(color = FALSE)
    }

    p
}



# Methods ======================================================================
#' @rdname plotGene
#' @export
setMethod(
    "plotGene",
    signature("matrix"),
    .plotGene)



#' @rdname plotGene
#' @export
setMethod(
    "plotGene",
    signature("bcbioRNASeq"),
    function(
        object,
        genes,
        normalized = "rlog",
        interestingGroups,
        medianLine = TRUE,
        color = scale_color_viridis(discrete = TRUE),
        return = "grid",
        headerLevel = 2L) {
        assert_is_a_string(normalized)
        if (missing(interestingGroups)) {
            interestingGroups <- bcbioBase::interestingGroups(object)
        }
        plotGene(
            object = counts(object, normalized = normalized),
            genes = genes,
            gene2symbol = gene2symbol(object),
            colData = colData(object),
            interestingGroups = interestingGroups,
            countsAxisLabel = paste(normalized, "counts"),
            medianLine = medianLine,
            color = color,
            return = return,
            headerLevel = headerLevel)
    })



#' @rdname plotGene
#' @export
setMethod(
    "plotGene",
    signature("DESeqDataSet"),
    function(
        object,
        genes,
        gene2symbol = NULL,
        interestingGroups = "sampleName",
        medianLine = TRUE,
        color = scale_color_viridis(discrete = TRUE),
        return = "grid",
        headerLevel = 2L) {
        plotGene(
            object = log2(counts(object, normalized = TRUE) + 1L),
            genes = genes,
            gene2symbol = gene2symbol,
            colData = colData(object),
            interestingGroups = interestingGroups,
            countsAxisLabel = "log2 counts",
            medianLine = medianLine,
            color = color,
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
        medianLine = TRUE,
        color = scale_color_viridis(discrete = TRUE),
        return = "grid",
        headerLevel = 2L) {
        if ("rlogIntercept" %in% colnames(mcols(object))) {
            normalized <- "rlog"
        } else {
            normalized <- "vst"
        }
        plotGene(
            object = assay(object),
            genes = genes,
            gene2symbol = gene2symbol,
            colData = colData(object),
            interestingGroups = interestingGroups,
            countsAxisLabel = paste(normalized, "counts"),
            medianLine = medianLine,
            color = color,
            return = return,
            headerLevel = headerLevel)
    })
