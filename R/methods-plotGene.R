#' Plot Individual Genes
#'
#' @name plotGene
#' @family Gene Expression Functions
#' @author Michael Steinbaugh
#'
#' @importFrom bcbioBase plotGene
#'
#' @inheritParams general
#' @param medianLine Include median line for each group. Disabled if samples
#'   are colored by sample name.
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
.plotGeneList <- function(
    object,
    countsAxisLabel = "counts",
    medianLine = TRUE,
    color = scale_color_viridis(discrete = TRUE)
) {
    stopifnot(is(object, "SummarizedExperiment"))
    stopifnot(length(rownames(object)) <= 50L)

    object <- convertGenesToSymbols(object)
    sampleData <- sampleData(object)
    interestingGroups <- interestingGroups(object)

    data <- .meltCounts(
        counts = assay(object),
        sampleData = sampleData(bcb_small)
    )

    # TODO Switch to DelayedArray approach here instead?
    list <- mclapply(
        X = rownames(object),
        FUN = function(x) {
            data <- data[data[["geneID"]] == x, , drop = FALSE]
            p <- ggplot(
                data = data,
                mapping = aes_string(
                    x = "geneID",
                    y = "counts",
                    color = "interestingGroups"
                )
            ) +
                genePoint() +
                theme(
                    axis.text.x = element_text(
                        angle = 90L, hjust = 1L, vjust = 0.5
                    )
                ) +
                labs(
                    title = geneID,
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
        }
    )
    names(list) <- rownames(object)
    list
}



.plotGeneFacet <- function(
    object,
    countsAxisLabel = "counts",
    medianLine = TRUE,
    color = scale_color_viridis(discrete = TRUE)
) {
    stopifnot(is(object, "SummarizedExperiment"))
    stopifnot(length(rownames(object)) <= 50L)

    object <- convertGenesToSymbols(object)
    sampleData <- sampleData(object)
    interestingGroups <- interestingGroups(object)

    data <- .meltCounts(
        counts = assay(object),
        sampleData = sampleData(bcb_small)
    )

    p <- ggplot(
        data = data,
        mapping = aes_string(
            x = "interestingGroups",
            y = "counts",
            color = "interestingGroups"
        )
    ) +
        genePoint() +
        theme(
            axis.text.x = element_text(angle = 90L, hjust = 1L, vjust = 0.5)
        ) +
        facet_wrap(facets = "geneID", scales = "free") +
        labs(
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



.plotGeneWide <- function(
    object,
    countsAxisLabel = "counts",
    medianLine = TRUE,
    color = scale_color_viridis(discrete = TRUE)
) {
    stopifnot(is(object, "SummarizedExperiment"))
    stopifnot(length(rownames(object)) <= 50L)

    object <- convertGenesToSymbols(object)
    sampleData <- sampleData(object)
    interestingGroups <- interestingGroups(object)

    data <- .meltCounts(
        counts = assay(object),
        sampleData = sampleData(bcb_small)
    )

    p <- ggplot(
        data = data,
        mapping = aes_string(
            x = "geneID",
            y = "counts",
            color = "interestingGroups"
        )
    ) +
        genePoint() +
        theme(
            axis.text.x = element_text(angle = 90L, hjust = 1L, vjust = 0.5)
        ) +
        labs(
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
#' @export
setMethod(
    "plotGene",
    signature("SummarizedExperiment"),
    function(
        object,
        genes,
        interestingGroups,
        countsAxisLabel = "counts",
        medianLine = TRUE,
        color = scale_color_viridis(discrete = TRUE),
        headerLevel = 2L,
        return = c("facet", "wide", "grid", "markdown", "list")
    ) {
        validObject(object)
        assert_is_a_bool(medianLine)
        assertIsColorScaleDiscreteOrNULL(color)
        assertIsAHeaderLevel(headerLevel)
        return <- match.arg(return)

        if (!missing(interestingGroups)) {
            interestingGroups(object) <- interestingGroups
        }

        rse <- as(object, "RangedSummarizedExperiment")
        rse <- rse[genes, , drop = FALSE]

        # Obtain ggplot objects per gene
        if (return %in% c("grid", "list", "markdown")) {
            plotlist <- .plotGeneList(
                object = rse,
                countsAxisLabel = countsAxisLabel,
                medianLine = medianLine,
                color = color
            )
        }

        if (return == "facet") {
            .plotGeneFacet(
                object = rse,
                countsAxisLabel = countsAxisLabel,
                medianLine = medianLine,
                color = color
            )
        } else if (return == "wide") {
            .plotGeneWide(
                object = rse,
                countsAxisLabel = countsAxisLabel,
                medianLine = medianLine,
                color = color
            )
        } else if (return == "grid") {
            if (length(plotlist) > 1L) {
                labels <- "AUTO"
            } else {
                labels <- NULL
            }
            plot_grid(plotlist = plotlist, labels = labels)
        } else if (return == "markdown") {
            markdownPlotlist(plotlist, headerLevel = headerLevel)
        } else if (return == "list") {
            plotlist
        }
    }
)



#' @rdname plotGene
#' @export
setMethod(
    "plotGene",
    signature("bcbioRNASeq"),
    function(
        object,
        normalized = c("rlog", "vst", "tpm"),
        ...
    ) {
        validObject(object)
        normalized <- match.arg(normalized)

        counts <- counts(object, normalized = normalized)
        # Ensure counts are log2 scale
        if (!normalized %in% c("rlog", "vst")) {
            counts <- log2(counts + 1L)
        }
        countsAxisLabel <- paste(normalized, "counts (log2)")

        # Coerce to RangedSummarizedExperiment and subset the genes
        rse <- as(object, "RangedSummarizedExperiment")
        assay(rse) <- counts

        # RangedSummarizedExperiment
        plotGene(
            object = rse,
            countsAxisLabel = countsAxisLabel,
            ...
        )
    }
)



#' @rdname plotGene
#' @export
setMethod(
    "plotGene",
    signature("DESeqDataSet"),
    function(object, ...) {
        validObject(object)

        # Ensure counts are log2 scale
        counts <- log2(counts(object, normalized = TRUE) + 1L)
        rse <- as(object, "RangedSummarizedExperiment")
        assay(rse) <- counts

        # RangedSummarizedExperiment
        plotGene(
            object = rse,
            countsAxisLabel = "normalized counts (log2)",
            ...
        )
    }
)



#' @rdname plotGene
#' @export
setMethod(
    "plotGene",
    signature("DESeqTransform"),
    function(object, ...) {
        validObject(object)
        if ("rlogIntercept" %in% colnames(mcols(object))) {
            normalized <- "rlog"
        } else {
            normalized <- "vst"
        }
        countsAxisLabel <- paste(normalized, "counts (log2)")
        rse <- as(object, "RangedSummarizedExperiment")
        plotGene(
            object = rse,
            countsAxisLabel = countsAxisLabel,
            ...
        )
    }
)
