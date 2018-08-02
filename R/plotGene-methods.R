#' Plot Individual Genes
#'
#' @name plotGene
#' @family Gene Expression Functions
#' @author Michael Steinbaugh
#'
#' @importFrom bcbioBase plotGene
#' @export
#'
#' @inheritParams general
#' @param countsAxisLabel `string`. Label to use for the counts axis.
#' @param medianLine `boolean`. Include median line for each group. Disabled if
#'   samples are colored by sample name.
#'
#' @return
#' - "`facet`": `ggplot` grouped by `sampleName`, with [ggplot2::facet_wrap()]
#'   applied to panel the samples.
#' - "`wide`": Show `ggplot` in wide format, with genes on the x-axis.
#'
#' @seealso [DESeq2::plotCounts()].
#'
#' @examples
#' # Gene identifiers
#' genes <- head(rownames(bcb_small), 4L)
#' glimpse(genes)
#'
#' # bcbioRNASeq ====
#' plotGene(bcb_small, genes = genes, return = "facet")
#' plotGene(bcb_small, genes = genes, return = "wide")
#'
#' # DESeqDataSet ====
#' plotGene(dds_small, genes = genes)
#'
#' # DESeqTransform ====
#' plotGene(rld_small, genes = genes)
NULL



# Constructors =================================================================
.plotGeneFacet <- function(
    object,
    countsAxisLabel = "counts",
    medianLine = TRUE,
    color = NULL,
    legend = TRUE
) {
    stopifnot(is(object, "SummarizedExperiment"))
    stopifnot(length(rownames(object)) <= 50L)

    object <- convertGenesToSymbols(object)
    interestingGroups <- .prepareInterestingGroups(
        object = object,
        interestingGroups = interestingGroups
    )

    data <- .meltCounts(
        counts = assay(object),
        sampleData = sampleData(object, clean = FALSE)
    )

    p <- ggplot(
        data = data,
        mapping = aes(
            x = !!sym("interestingGroups"),
            y = !!sym("counts"),
            color = !!sym("interestingGroups")
        )
    ) +
        .genePoint(show.legend = legend) +
        facet_wrap(facets = sym("geneID"), scales = "free_y") +
        labs(
            x = NULL,
            y = countsAxisLabel,
            color = paste(interestingGroups, collapse = ":\n")
        )

    if (isTRUE(medianLine) && !identical(interestingGroups, "sampleName")) {
        p <- p + .geneMedianLine
    }

    if (is(color, "ScaleDiscrete")) {
        p <- p + color
    }

    if (identical(interestingGroups, "sampleName")) {
        p <- p + guides(color = FALSE)
    }

    p
}



.plotGeneList <- function(
    object,
    countsAxisLabel = "counts",
    medianLine = TRUE,
    color = NULL,
    legend = TRUE
) {
    stopifnot(is(object, "SummarizedExperiment"))
    stopifnot(length(rownames(object)) <= 50L)

    object <- convertGenesToSymbols(object)
    interestingGroups <- .prepareInterestingGroups(
        object = object,
        interestingGroups = interestingGroups
    )

    data <- .meltCounts(
        counts = assay(object),
        sampleData = sampleData(object, clean = FALSE)
    )

    list <- lapply(
        X = rownames(object),
        FUN = function(geneID) {
            data <- data[data[["geneID"]] == geneID, , drop = FALSE]
            p <- ggplot(
                data = data,
                mapping = aes(
                    x = !!sym("geneID"),
                    y = !!sym("counts"),
                    color = !!sym("interestingGroups")
                )
            ) +
                .genePoint(show.legend = legend) +
                labs(
                    title = geneID,
                    x = NULL,
                    y = countsAxisLabel,
                    color = paste(interestingGroups, collapse = ":\n")
                ) +
                theme(
                    axis.text.x = element_blank(),
                    axis.ticks.x = element_blank()
                )

            if (
                isTRUE(medianLine) &&
                !identical(interestingGroups, "sampleName")
            ) {
                p <- p + .geneMedianLine
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



.plotGeneWide <- function(
    object,
    countsAxisLabel = "counts",
    medianLine = TRUE,
    color = NULL,
    legend = TRUE
) {
    stopifnot(is(object, "SummarizedExperiment"))
    stopifnot(length(rownames(object)) <= 50L)

    object <- convertGenesToSymbols(object)
    interestingGroups <- .prepareInterestingGroups(
        object = object,
        interestingGroups = interestingGroups
    )

    data <- .meltCounts(
        counts = assay(object),
        sampleData = sampleData(object, clean = FALSE)
    )

    p <- ggplot(
        data = data,
        mapping = aes(
            x = !!sym("geneID"),
            y = !!sym("counts"),
            color = !!sym("interestingGroups")
        )
    ) +
        .genePoint(show.legend = legend) +
        labs(
            x = NULL,
            y = countsAxisLabel,
            color = paste(interestingGroups, collapse = ":\n")
        )

    if (isTRUE(medianLine) && !identical(interestingGroups, "sampleName")) {
        p <- p + .geneMedianLine
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
        color = getOption("bcbio.discrete.color", NULL),
        legend = getOption("bcbio.legend", TRUE),
        headerLevel = 2L,
        return = c("facet", "wide")
    ) {
        validObject(object)
        assert_is_character(genes)
        interestingGroups <- .prepareInterestingGroups(
            object = object,
            interestingGroups = interestingGroups
        )
        assert_is_a_bool(medianLine)
        assertIsColorScaleDiscreteOrNULL(color)
        assert_is_a_bool(legend)
        assertIsAHeaderLevel(headerLevel)
        return <- match.arg(return)

        # Coerce to RangedSummarizedExperiment
        rse <- object %>%
            as("RangedSummarizedExperiment") %>%
            .[genes, , drop = FALSE]

        if (return == "facet") {
            fxn <- .plotGeneFacet
        } else if (return == "wide") {
            fxn <- .plotGeneWide
        }

        fxn(
            object = rse,
            countsAxisLabel = countsAxisLabel,
            medianLine = medianLine,
            color = color,
            legend = legend
        )
    }
)



#' @rdname plotGene
#' @export
setMethod(
    "plotGene",
    signature("bcbioRNASeq"),
    function(
        object,
        normalized = c("vst", "rlog", "tmm", "tpm", "rle"),
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
