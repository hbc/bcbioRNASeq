#' Sample PCA plot for transformed data
#'
#' Wrapper for [DESeq2::plotPCA()] that improves principal component analysis
#' (PCA) sample coloring and labeling.
#'
#' PCA (Jolliffe, et al., 2002) is a multivariate technique that allows us to
#' summarize the systematic patterns of variations in the data. PCA takes the
#' expression levels for genes and transforms it in principal component space,
#' reducing each sample into one point. Thereby, we can separate samples by
#' expression variation, and identify potential sample outliers. The PCA plot is
#' a way to look at how samples are clustering.
#'
#' @name plotPCA
#' @family Quality Control Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @param ntop `scalar integer` or `Inf`. Number of most variable genes to plot.
#'   Use `Inf` to include all genes.
#' @param ... Additional arguments.
#'
#' @seealso
#' - [DESeq2::plotPCA()].
#' - `getMethod("plotPCA", "DESeqTransform")`
#'
#' @return `ggplot` or `data.frame`.
#'
#' @examples
#' # bcbioRNASeq ====
#' plotPCA(
#'     object = bcb_small,
#'     normalized = "vst",
#'     label = TRUE
#' )
#' plotPCA(
#'     object = bcb_small,
#'     normalized = "rlog",
#'     interestingGroups = "sampleName",
#'     label = FALSE
#' )
NULL



#' @rdname plotPCA
#' @name plotPCA
#' @importFrom BiocGenerics plotPCA
#' @usage plotPCA(object, ...)
#' @export
NULL



plotPCA.SummarizedExperiment <-  # nolint
    function(
        object,
        interestingGroups,
        ntop = 500L,
        color = getOption("bcbio.discrete.color", NULL),
        label = getOption("bcbio.label", FALSE),
        title = "pca",
        subtitle = NULL,
        return = c("ggplot", "data.frame")
    ) {
        # Legacy arguments -----------------------------------------------------
        # nocov start
        call <- match.call()
        # genes
        if ("genes" %in% names(call)) {
            stop("`genes` is defunct. Use `ntop` argument instead.")
        }
        # samples, censorSamples
        if (any(c("samples", "censorSamples") %in% names(call))) {
            stop("Sample selection is defunct. Use bracket-based subsetting.")
        }
        # returnData
        if ("returnData" %in% names(call)) {
            warning("`returnData` is deprecated in favor of `return`")
            returnData <- call[["returnData"]]
            if (isTRUE(returnData)) {
                return <- "data.frame"
            }
        }
        # nocov end

        # Assert checks --------------------------------------------------------
        validObject(object)
        interestingGroups <- matchInterestingGroups(
            object = object,
            interestingGroups = interestingGroups
        )
        assert_is_a_number(ntop)
        assertIsColorScaleDiscreteOrNULL(color)
        assert_is_a_bool(label)
        assertIsAStringOrNULL(title)
        # exists check is needed for legacy `returnData` match above
        if (!exists(return, inherits = FALSE)) {
            return <- match.arg(return)
        }

        if (identical(ntop, Inf)) {
            nGene <- nrow(object)
        } else {
            nGene <- ntop
        }
        message(paste("Plotting PCA using", nGene, "genes"))

        # Get coordinates using DESeqTransform method
        dt <- DESeqTransform(object)
        colData <- colData(object)
        data <- plotPCA(
            object = dt,
            intgroup = interestingGroups,
            ntop = ntop,
            returnData = TRUE
        ) %>%
            camel()

        # Early return if `data.frame` is requested
        if (return == "data.frame") {
            return(data)
        }

        # Use `sampleName` for plot labels
        if (isTRUE(label)) {
            data[["sampleName"]] <- colData[["sampleName"]]
        }

        percentVar <- round(100L * attr(data, "percentVar"))

        # `DESeq2::plotPCA()` defines interesting groups in `group` column
        p <- ggplot(
            data = data,
            mapping = aes(
                x = !!sym("pc1"),
                y = !!sym("pc2"),
                color = !!sym("group")
            )
        ) +
            geom_point(size = 4L) +
            coord_fixed() +
            labs(
                title = title,
                subtitle = subtitle,
                x = paste0("pc1: ", percentVar[[1L]], "% variance"),
                y = paste0("pc2: ", percentVar[[2L]], "% variance"),
                color = paste(interestingGroups, collapse = ":\n")
            )

        if (is(color, "ScaleDiscrete")) {
            p <- p + color
        }

        if (isTRUE(label)) {
            p <- p + bcbio_geom_label_repel(
                mapping = aes(label = !!sym("sampleName"))
            )
        }

        p
    }



#' @rdname plotPCA
#' @export
setMethod(
    f = "plotPCA",
    signature = signature("SummarizedExperiment"),
    definition = plotPCA.SummarizedExperiment
)



plotPCA.bcbioRNASeq <-  # nolint
    function(
        object,
        normalized = c("vst", "rlog", "tmm", "tpm", "rle"),
        ...
    ) {
        validObject(object)
        normalized <- match.arg(normalized)
        message(paste("Using", normalized, "counts"))
        counts <- counts(object, normalized = normalized)
        rse <- as(object, "RangedSummarizedExperiment")
        assay(rse) <- counts
        plotPCA(rse, ...)
    }



#' @rdname plotPCA
#' @export
setMethod(
    f = "plotPCA",
    signature = signature("bcbioRNASeq"),
    definition = plotPCA.bcbioRNASeq
)
