#' Extract or replace parts of an object
#'
#' Extract genes by row and samples by column from a `bcbioRNASeq` object.
#' Internal count transformations are rescaled automatically, if defined.
#'
#' @note DESeq2 transformations will only be calculated when `transform = TRUE`
#'   and either `rlog` or `vst` counts are defined in [assays()].
#'
#' @name extract
#' @family S4 Object
#' @author Lorena Pantano, Michael Steinbaugh
#'
#' @inheritParams base::`[`
#' @inheritParams general
#'
#' @param transform `boolean`. Recalculate slotted DESeq2 transformations.
#'   Recommended by default but can take a while for large datasets.
#'
#' @return `bcbioRNASeq`.
#'
#' @seealso `help("[", "base")`.
#'
#' @examples
#' object <- bcb_small
#'
#' # Minimum of 100 genes, 2 samples
#' genes <- head(rownames(object), 100L)
#' head(genes)
#' samples <- head(colnames(object), 2L)
#' head(samples)
#'
#' # Extract by sample name
#' object[, samples]
#'
#' # Extract by gene list
#' object[genes, ]
#'
#' # Extract by both genes and samples
#' x <- object[genes, samples]
#' print(x)
#' assayNames(x)
#'
#' # Fast subsetting, by skipping DESeq2 transformations
#' # Note that rlog and vst counts will not be defined
#' x <- object[, samples, transform = FALSE]
#' print(x)
#' names(assays(x))
NULL



#' @rdname extract
#' @export
setMethod(
    f = "[",
    signature = signature(
        x = "bcbioRNASeq",
        i = "ANY",
        j = "ANY",
        drop = "ANY"
    ),
    definition = function(
        x,
        i,
        j,
        drop = FALSE,
        transform = TRUE
    ) {
        validObject(x)

        # Genes (rows)
        if (missing(i)) {
            i <- 1L:nrow(x)
        }
        # Require at least 100 genes
        assert_all_are_in_range(length(i), lower = 100L, upper = Inf)

        # Samples (columns)
        if (missing(j)) {
            j <- 1L:ncol(x)
        }
        # Require at least 2 samples
        assert_all_are_in_range(length(j), lower = 2L, upper = Inf)

        # Regenerate RangedSummarizedExperiment
        rse <- as(x, "RangedSummarizedExperiment")
        rse <- rse[i, j, drop = drop]

        # Early return if dimensions are unmodified
        if (identical(dim(rse), dim(x))) {
            return(x)
        }

        # Assays ---------------------------------------------------------------
        assays <- assays(rse)

        # Always recalculate DESeq2 normalized counts
        message("Updating normalized counts")
        dds <- .regenerateDESeqDataSet(rse)
        assays[["normalized"]] <- counts(dds, normalized = TRUE)

        # Optionally, recalculate DESeq2 transformations (rlog, vst)
        if (
            any(c("rlog", "vst") %in% names(assays)) &&
            isTRUE(transform)
        ) {
            # Update DESeq2 transformations by default, but only if they are
            # already defined in assays (rlog, vst)
            message(paste(
                "Recalculating DESeq2 variance stabilizations",
                "(transform = TRUE)"
            ))
            # vst
            if ("vst" %in% names(assays)) {
                message("Applying variance-stabilizing transformation")
                assays[["vst"]] <- assay(varianceStabilizingTransformation(dds))
            }
            # rlog
            if ("rlog" %in% names(assays)) {
                message("Applying rlog transformation")
                assays[["rlog"]] <- assay(rlog(dds))
            }
        } else {
            # Otherwise, ensure previous calculations are removed from assays
            message(paste(
                "Skipping DESeq2 variance stabilizations",
                "(transform = FALSE)"
            ))
            assays[["rlog"]] <- NULL
            assays[["vst"]] <- NULL
        }

        # Column data ----------------------------------------------------------
        # Ensure factors get releveled
        colData <- colData(rse) %>%
            as.data.frame() %>%
            rownames_to_column() %>%
            mutate_if(is.character, as.factor) %>%
            mutate_if(is.factor, droplevels) %>%
            column_to_rownames() %>%
            as("DataFrame")

        # Metadata -------------------------------------------------------------
        metadata <- metadata(rse)
        metadata[["subset"]] <- TRUE
        # Update version, if necessary
        if (!identical(metadata[["version"]], packageVersion)) {
            metadata[["originalVersion"]] <- metadata[["version"]]
            metadata[["version"]] <- packageVersion
        }

        # Return ---------------------------------------------------------------
        .new.bcbioRNASeq(
            assays = assays,
            rowRanges = rowRanges(rse),
            colData = colData,
            metadata = metadata
        )
    }
)
