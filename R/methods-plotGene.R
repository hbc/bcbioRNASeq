#' Plot Individual Genes
#'
#' @name plotGene
#' @family Gene Expression Functions
#' @author Michael Steinbaugh
#'
#' @importFrom bcbioBase plotGene
#'
#' @inheritParams general
#' @param medianLine Include median line for each group. Disabled when samples
#'   are grouped by `sampleName`.
#'
#' @return
#' - "`grid`": Show [cowplot::plot_grid()], paneled per gene.
#' - "`wide`": Show `ggplot` in wide format, with genes on the x-axis.
#' - "`list`": `list`, containing per gene `ggplot` objects.
#' - "`markdown`": Show tabset R Markdown output, tabbed per gene.
#'
#' @seealso [DESeq2::plotCounts()].
#'
#' @examples
#' load(system.file("extdata/bcb_small.rda", package = "bcbioRNASeq"))
#' load(system.file("extdata/dds_small.rda", package = "bcbioRNASeq"))
#' load(system.file("extdata/rld_small.rda", package = "bcbioRNASeq"))
#'
#' # Gene identifiers
#' genes <- head(rownames(bcb_small), 8L)
#'
#' # bcbioRNASeq ====
#' plotGene(bcb_small, genes = genes, return = "grid")
#' plotGene(bcb_small, genes = genes, return = "wide")
#'
#' # DESeqDataSet ====
#' plotGene(dds_small, genes = genes)
#'
#' # DESeqTransform ====
#' plotGene(rld_small, genes = genes)
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
    headerLevel = 2L,
    return = c("grid", "wide", "list", "markdown")
) {
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
    assertIsAHeaderLevel(headerLevel)
    return <- match.arg(return)

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
            color = color
        )
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
            color = color
        )
    } else if (return == "list") {
        plotlist
    } else if (return == "markdown") {
        markdownPlotlist(plotlist, headerLevel = headerLevel)
    }
}



#' @importFrom ggplot2 aes_string element_text expand_limits geom_point ggplot
#'   guides labs theme
#' @importFrom parallel mcmapply
#' @importFrom tibble tibble
.plotGeneList <- function(
    object,
    genes,
    gene2symbol = NULL,
    colData,
    interestingGroups = "sampleName",
    countsAxisLabel = "log2 counts",
    medianLine = TRUE,
    color = scale_color_viridis(discrete = TRUE)
) {
    assert_is_subset(genes, rownames(object))
    object <- object[genes, , drop = FALSE]

    # Gene to symbol mappings
    if (is.data.frame(gene2symbol)) {
        assert_is_subset(genes, gene2symbol[["geneID"]])
        match <- match(x = genes, table = gene2symbol[["geneID"]])
        symbols <- gene2symbol[match, "geneName", drop = TRUE]
        rownames(object) <- symbols
        genes <- symbols
    }

    mcmapply(
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
                interestingGroups = colData[["interestingGroups"]]
            )

            p <- ggplot(
                data = data,
                mapping = aes_string(
                    x = "x",
                    y = "y",
                    color = "interestingGroups"
                )
            ) +
                genePoint() +
                theme(axis.text.x = element_text(angle = 90L)) +
                labs(
                    title = gene,
                    x = NULL,
                    y = countsAxisLabel,
                    color = paste(interestingGroups, collapse = ":\n")
                ) +
                # expand_limits(y = 0L) +
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
        USE.NAMES = TRUE
    )
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
    countsAxisLabel = "log2 counts",
    medianLine = TRUE,
    color = scale_color_viridis(discrete = TRUE),
    title = NULL
) {
    # Gene to symbol mappings
    if (is.data.frame(gene2symbol)) {
        assert_is_subset(genes, gene2symbol[["geneID"]])
        match <- match(x = genes, table = gene2symbol[["geneID"]])
        symbols <- gene2symbol[match, "geneName", drop = TRUE]
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
        left_join(as.data.frame(colData), by = "sampleID")

    p <- ggplot(
        data = data,
        mapping = aes_string(
            x = "gene",
            y = "counts",
            color = "interestingGroups"
        )
    ) +
        genePoint() +
        theme(axis.text.x = element_text(angle = 90L, hjust = 1L)) +
        labs(
            title = title,
            x = NULL,
            y = countsAxisLabel,
            color = paste(interestingGroups, collapse = ":\n")
        )

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
#' @importFrom bcbioBase gene2symbol
#' @export
setMethod(
    "plotGene",
    signature("bcbioRNASeq"),
    function(
        object,
        genes,
        normalized = c("rlog", "vst", "tpm"),
        interestingGroups,
        medianLine = TRUE,
        color = scale_color_viridis(discrete = TRUE),
        headerLevel = 2L,
        return = c("grid", "wide", "list", "markdown")
    ) {
        # Passthrough: genes, medianLine, color, headerLevel, return
        normalized <- match.arg(normalized)
        if (missing(interestingGroups)) {
            interestingGroups <- bcbioBase::interestingGroups(object)
        }
        counts <- counts(object, normalized = normalized)
        # Ensure counts are log2 scale
        if (!normalized %in% c("rlog", "vst")) {
            counts <- log2(counts + 1L)
        }
        .plotGene(
            object = counts,
            genes = genes,
            gene2symbol = gene2symbol(object),
            colData = colData(object),
            interestingGroups = interestingGroups,
            countsAxisLabel = paste(normalized, "counts (log2)"),
            medianLine = medianLine,
            color = color,
            return = return,
            headerLevel = headerLevel
        )
    }
)



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
        headerLevel = 2L,
        return = c("grid", "wide", "list", "markdown")
    ) {
        # Passthrough: genes, gene2symbol, medianLine, color, headerLevel,
        # return
        counts <- log2(counts(object, normalized = TRUE) + 1L)
        colData <- colData(object)
        colData[["sizeFactor"]] <- NULL
        .plotGene(
            object = counts,
            genes = genes,
            gene2symbol = gene2symbol,
            colData = colData,
            interestingGroups = interestingGroups,
            countsAxisLabel = "normalized counts (log2)",
            medianLine = medianLine,
            color = color,
            return = return,
            headerLevel = headerLevel
        )
    }
)



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
        headerLevel = 2L,
        return = c("grid", "wide", "list", "markdown")
    ) {
        # Passthrough: genes, gene2symbol, medianLine, color, headerLevel,
        # return
        if ("rlogIntercept" %in% colnames(mcols(object))) {
            normalized <- "rlog"
        } else {
            normalized <- "vst"
        }
        counts <- assay(object)
        colData <- colData(object)
        colData[["sizeFactor"]] <- NULL
        .plotGene(
            object = counts,
            genes = genes,
            gene2symbol = gene2symbol,
            colData = colData,
            interestingGroups = interestingGroups,
            countsAxisLabel = paste(normalized, "counts (log2)"),
            medianLine = medianLine,
            color = color,
            return = return,
            headerLevel = headerLevel
        )
    }
)
