#' Extract or Replace Parts of an Object
#'
#' Extract genes by row and samples by column from a `bcbioRNASeq` object.
#' Internal count transformations are rescaled automatically, if defined.
#'
#' @note DESeq2 transformations will only be updated when `recalculate = TRUE`
#'   and either `rlog` or `vst` counts are defined in [assays()].
#'
#' @name extract
#' @family S4 Functions
#' @author Lorena Pantano, Michael Steinbaugh
#'
#' @inheritParams base::`[`
#' @inheritParams general
#'
#' @param recalculate `boolean`. Recalculate DESeq2 normalized counts and
#'   variance-stabilizing transformations defined in [assays()]. Recommended by
#'   default, but can take a long time for large datasets.
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
#' # Fast subsetting, by skipping DESeq2 recalculations.
#' # Note that `normalized`, `rlog`, and `vst` assays will not be defined.
#' x <- object[, samples, recalculate = FALSE]
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
        drop = "ANY"  # Don't use logical here.
    ),
    definition = function(
        x,
        i,
        j,
        drop = FALSE,
        recalculate = TRUE
    ) {
        validObject(x)
        # Never allow the user to drop on extraction.
        stopifnot(!isTRUE(drop))
        assert_is_a_bool(recalculate)

        # Genes (rows)
        if (missing(i)) {
            i <- seq_len(nrow(x))
        }
        # Require at least 100 genes.
        assert_all_are_in_range(length(i), lower = 100L, upper = Inf)

        # Samples (columns)
        if (missing(j)) {
            j <- seq_len(ncol(x))
        }
        # Require at least 2 samples.
        assert_all_are_in_range(length(j), lower = 2L, upper = Inf)

        # Early return if dimensions are unmodified.
        if (identical(
            x = dim(x),
            y = c(length(i), length(j))
        )) {
            message("Returning object unmodified.")
            return(x)
        }

        # Regenerate RangedSummarizedExperiment.
        rse <- as(x, "RangedSummarizedExperiment")
        rse <- rse[i, j, drop = FALSE]

        # Metadata -------------------------------------------------------------
        metadata <- metadata(rse)
        metadata[["subset"]] <- TRUE

        # Resize inferential replicates, if defined.
        infReps <- metadata[["infReps"]]
        if (
            is.list(infReps) &&
            has_length(infReps)
        ) {
            message("Resizing inferential replicates...")
            infReps <- infReps[colnames(rse)]
            infReps <- lapply(infReps, function(x) {
                x[rownames(rse), , drop = FALSE]
            })
        } else {
            infReps <- NULL
        }
        metadata[["infReps"]] <- infReps
        metadata(rse) <- metadata

        # Assays ---------------------------------------------------------------
        assays <- assays(rse)

        # Recalculate DESeq2 normalized counts and variance stabilizations.
        if (isTRUE(recalculate)) {
            message(paste(
                "Recalculating DESeq2 normalizations",
                "(recalculate = TRUE)..."
            ))

            # Normalized counts
            dds <- .regenerateDESeqDataSet(rse)
            assays[["normalized"]] <- counts(dds, normalized = TRUE)

            # Variance-stabilizing transformations.
            if (any(c("rlog", "vst") %in% names(assays))) {
                # Update DESeq2 transformations by default, but only if they are
                # already defined in assays (rlog, vst).
                message("Recalculating DESeq2 variance stabilizations...")
                # vst
                if ("vst" %in% names(assays)) {
                    message("Applying variance-stabilizing transformation...")
                    assays[["vst"]] <-
                        assay(varianceStabilizingTransformation(dds))
                }
                # rlog
                if ("rlog" %in% names(assays)) {
                    message("Applying rlog transformation...")
                    assays[["rlog"]] <- assay(rlog(dds))
                }
            }
        } else {
            # Otherwise, ensure previous calculations are removed from assays.
            message(paste(
                "Skipping DESeq2 normalizations",
                "(recalculate = FALSE)"
            ))
            assays[["normalized"]] <- NULL
            assays[["rlog"]] <- NULL
            assays[["vst"]] <- NULL
        }

        # Row data -------------------------------------------------------------
        # Ensure factors get releveled.
        rowRanges <- rowRanges(rse)
        mcols(rowRanges) <- mcols(rowRanges) %>%
            as("tbl_df") %>%
            mutate_if(is.factor, droplevels) %>%
            as("DataFrame")

        # Column data ----------------------------------------------------------
        # Ensure factors get releveled.
        colData <- colData(rse) %>%
            as.data.frame() %>%
            rownames_to_column() %>%
            mutate_if(is.character, as.factor) %>%
            mutate_if(is.factor, droplevels) %>%
            column_to_rownames() %>%
            as("DataFrame")

        # Return ---------------------------------------------------------------
        .new.bcbioRNASeq(
            assays = assays,
            rowRanges = rowRanges,
            colData = colData,
            metadata = metadata
        )
    }
)
