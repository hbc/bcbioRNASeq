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
#'
#' @importFrom BiocGenerics plotPCA
#' @export
#'
#' @inheritParams general
#'
#' @seealso
#' - [DESeq2::plotPCA()].
#' - `getMethod("plotPCA", "DESeqTransform")`
#'
#' @return `ggplot` or `data.frame`.
#'
#' @examples
#' # bcbioRNASeq ====
#' plotPCA(bcb_small, label = FALSE)
#' plotPCA(bcb_small, label = TRUE)
#'
#' # Select samples
#' plotPCA(
#'     object = bcb_small,
#'     samples = head(colnames(bcb_small), 4L),
#'     label = TRUE
#' )
#'
#' # DESeqDataSet ====
#' # DESeqTransform method is preferred
#' plotPCA(dds_small)
NULL



# Methods ======================================================================
#' @rdname plotPCA
#' @export
setMethod(
    "plotPCA",
    signature("SummarizedExperiment"),
    function(
        object,
        genes = NULL,
        samples = NULL,
        interestingGroups,
        color = NULL,
        label = FALSE,
        title = "pca",
        subtitle = NULL,
        return = c("ggplot", "data.frame")
    ) {
        # Legacy arguments =====================================================
        call <- match.call()
        # censorSamples
        if ("censorSamples" %in% names(call)) {
            warning("`censorSamples` is deprecated in favor of `samples`")
            censorSamples <- call[["censorSamples"]]
            assert_is_subset(censorSamples, colnames(object))
            samples <- setdiff(colnames(object), censorSamples)
        }
        # returnData
        if ("returnData" %in% names(call)) {
            warning("`returnData` is deprecated in favor of `return`")
            returnData <- call[["returnData"]]
            if (isTRUE(returnData)) {
                return <- "data.frame"
            }
        }

        # Assert checks ========================================================
        validObject(object)
        assert_is_any_of(genes, c("character", "NULL"))
        assert_is_any_of(samples, c("character", "NULL"))
        if (missing(interestingGroups)) {
            interestingGroups <- basejump::interestingGroups(object)
        } else {
            interestingGroups(object) <- interestingGroups
        }
        assertIsColorScaleDiscreteOrNULL(color)
        assert_is_a_bool(label)
        assertIsAStringOrNULL(title)
        # exists check is needed for legacy `returnData` match above
        if (!exists(return, inherits = FALSE)) {
            return <- match.arg(return)
        }

        # Prepare counts matrix ================================================
        # Subset samples, if desired
        if (length(samples)) {
            object <- object[, samples]
        }

        # Subset genes, if desired
        if (length(genes)) {
            object <- object[genes, , drop = FALSE]
            # Set ntop to the number of genes requested
            ntop <- length(genes)
            message(paste("Plotting PCA using", ntop, "genes"))
        } else {
            # Recommended DESeq default of most variable genes
            ntop <- 500L
            message(paste("Plotting PCA using", ntop, "most variable genes"))
        }

        # Get PCA data using DESeqTransform method =============================
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
)



#' @rdname plotPCA
#' @export
setMethod(
    "plotPCA",
    signature("bcbioRNASeq"),
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
)
