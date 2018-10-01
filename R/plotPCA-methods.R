#' Sample PCA Plot for Transformed Data
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
#' @importFrom BiocGenerics plotPCA
#' @export
#'
#' @inheritParams general
#' @param ntop `scalar integer` or `Inf`. Number of most variable genes to plot.
#'   Use `Inf` to include all genes.
#'
#' @seealso
#' - [DESeq2::plotPCA()].
#' - `getMethod("plotPCA", "DESeqTransform")`
#'
#' @return `ggplot` or `DataFrame`.
#'
#' @examples
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



# We're using this internal constructor for `plotDEGHeatmap()` also.
# Keep this code separate from the bcbioRNASeq method, but don't export.
.plotPCA.SummarizedExperiment <-  # nolint
    function(
        object,
        interestingGroups = NULL,
        ntop = 500L,
        color = getOption("basejump.discrete.color", NULL),
        label = getOption("basejump.label", FALSE),
        title = "PCA",
        subtitle = NULL,
        return = c("ggplot", "DataFrame")
    ) {
        # Legacy arguments -----------------------------------------------------
        # nocov start
        call <- standardizeCall()
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
            warning("`returnData` is deprecated in favor of `return`.")
            returnData <- call[["returnData"]]
            if (isTRUE(returnData)) {
                return <- "DataFrame"
            }
        }
        # nocov end

        validObject(object)
        interestingGroups <- matchInterestingGroups(
            object = object,
            interestingGroups = interestingGroups
        )
        interestingGroups(object) <- interestingGroups
        assert_is_a_number(ntop)
        assertIsColorScaleDiscreteOrNULL(color)
        assert_is_a_bool(label)
        assertIsAStringOrNULL(title)
        # `exists()` check here is needed for legacy `returnData` match above.
        if (!exists(return, inherits = FALSE)) {
            return <- match.arg(return)
        }

        if (identical(ntop, Inf)) {
            nGene <- nrow(object)
        } else {
            nGene <- ntop
        }
        message(paste("Plotting PCA using", nGene, "genes."))

        # Using the `DESeq2::plotPCA()` `DESeqTransform` method to obtain the
        # PCA coordinates in a data frame.
        data <- do.call(
            what = plotPCA,
            args = list(
                object = DESeqTransform(object),
                intgroup = interestingGroups,
                ntop = ntop,
                returnData = TRUE
            )
        )
        # Check to make sure data is returning as expected.
        assert_is_subset(
            x = c("PC1", "PC2", "group", "name"),
            y = colnames(data)
        )
        assert_are_identical(
            x = rownames(data),
            y = rownames(colData(object))
        )

        # Standardize the columns to match our conventions.
        colnames(data) <- gsub("^group$", "interestingGroups", colnames(data))
        data[["name"]] <- NULL
        data[["sampleName"]] <- colData(object)[["sampleName"]]

        # Make sure `percentVar` attribute is set.
        assert_is_numeric(attr(data, "percentVar"))

        # Early return if `DataFrame` is requested.
        if (return == "DataFrame") {
            return(as(data, "DataFrame"))
        }

        # Improve the appearance of the percent variation.
        percentVar <- round(100L * attr(data, "percentVar"))

        # `DESeq2::plotPCA()` defines interesting groups in `group` column.
        p <- ggplot(
            data = data,
            mapping = aes(
                x = !!sym("PC1"),
                y = !!sym("PC2"),
                color = !!sym("interestingGroups")
            )
        ) +
            geom_point(size = 4L) +
            coord_fixed() +
            labs(
                title = title,
                subtitle = subtitle,
                x = paste0("PC1: ", percentVar[[1L]], "% variance"),
                y = paste0("PC2: ", percentVar[[2L]], "% variance"),
                color = paste(interestingGroups, collapse = ":\n")
            )

        if (is(color, "ScaleDiscrete")) {
            p <- p + color
        }

        if (isTRUE(label)) {
            p <- p + basejump_geom_label_repel(
                mapping = aes(label = !!sym("sampleName"))
            )
        }

        p
    }



.plotPCA.bcbioRNASeq <-  # nolint
    function(
        object,
        normalized = c("vst", "rlog", "tmm", "tpm", "rle")
    ) {
        validObject(object)
        normalized <- match.arg(normalized)
        message(paste("Using", normalized, "counts."))
        rse <- as(object, "RangedSummarizedExperiment")
        assays(rse) <- list(counts(object, normalized = normalized))
        do.call(
            what = .plotPCA.SummarizedExperiment,
            args = matchArgsToDoCall(
                args = list(object = rse),
                removeFormals = "normalized"
            )
        )
    }
f1 <- formals(.plotPCA.bcbioRNASeq)
f2 <- formals(.plotPCA.SummarizedExperiment)
f2 <- f2[setdiff(names(f2), names(f1))]
f <- c(f1, f2)
formals(.plotPCA.bcbioRNASeq) <- f



#' @rdname plotPCA
#' @export
setMethod(
    f = "plotPCA",
    signature = signature("bcbioRNASeq"),
    definition = .plotPCA.bcbioRNASeq
)
