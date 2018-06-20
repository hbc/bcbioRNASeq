#' Sample PCA Plot for Transformed Data
#'
#' Wrapper for [DESeq2::plotPCA()] that improves principal component analysis
#' (PCA) sample coloring and labeling.
#'
#' @name plotPCA
#' @family Quality Control Functions
#' @author Michael Steinbaugh
#'
#' @importFrom BiocGenerics plotPCA
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
        title = "PCA",
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
            interestingGroups <- bcbioBase::interestingGroups(object)
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
        colData <- colData(dt)
        assert_is_subset("sampleName", colnames(colData))
        data <- plotPCA(
            object = dt,
            intgroup = interestingGroups,
            ntop = ntop,
            returnData = TRUE
        )
        data[["sampleName"]] <- colData[["sampleName"]]

        # Early return if `data.frame` is requested
        if (return == "data.frame") {
            return(data)
        }

        percentVar <- round(100L * attr(data, "percentVar"))

        # `DESeq2::plotPCA()` defines interesting groups in `group` column
        p <- ggplot(
            data = data,
            mapping = aes_string(
                x = "PC1",
                y = "PC2",
                color = "group"
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
            p <- p + bcbio_geom_label_repel(
                mapping = aes_string(label = "sampleName")
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
