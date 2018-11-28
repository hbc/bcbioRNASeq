#' @name extract
#' @author Michael Steinbaugh, Lorena Pantano
#' @inherit base::Extract title params references
#'
#' @description Extract genes by row and samples by column.
#'
#' @details
#' Internal count transformations are rescaled automatically, if defined. DESeq2
#' transformations will only be updated when `recalculate = TRUE` and either
#' `rlog` or `vst` counts are defined in `assays()`.
#
#' @inheritParams params
#' @param recalculate `boolean`. Recalculate DESeq2 normalized counts and
#'   variance-stabilizing transformations defined in `assays()`. Recommended by
#'   default, but can take a long time for large datasets.
#'
#' @return `bcbioRNASeq`.
#'
#' @examples
#' data(bcb)
#' object <- bcb
#'
#' ## Minimum of 100 genes, 2 samples.
#' genes <- head(rownames(object), 100L)
#' head(genes)
#' samples <- head(colnames(object), 2L)
#' head(samples)
#'
#' ## Extract by sample name.
#' object[, samples]
#'
#' ## Extract by gene list.
#' object[genes, ]
#'
#' ## Extract by both genes and samples.
#' x <- object[genes, samples]
#' print(x)
#' assayNames(x)
#'
#' ## Fast subsetting, by skipping DESeq2 recalculations.
#' ## Note that `normalized`, `rlog`, and `vst` assays will not be defined.
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
        assert_that(!isTRUE(drop))
        assert_is_a_bool(recalculate)

        # Genes (rows)
        if (missing(i)) {
            i <- seq_len(nrow(x))
        }
        # Require at least 50 genes.
        assert_all_are_in_range(length(i), lower = 50L, upper = Inf)

        # Samples (columns)
        if (missing(j)) {
            j <- seq_len(ncol(x))
        }
        # Require at least 2 samples.
        assert_all_are_in_range(length(j), lower = 2L, upper = Inf)

        # Don't attempt to recalculate normalizations if the dimensions remain
        # unchanged.
        if (identical(
            x = dim(x),
            y = c(length(i), length(j))
        )) {
            subset <- FALSE
        } else {
            subset <- TRUE
        }

        # Regenerate RangedSummarizedExperiment.
        rse <- as(x, "RangedSummarizedExperiment")
        rse <- rse[i, j, drop = FALSE]

        # Assays ---------------------------------------------------------------
        assays <- assays(rse)

        if (isTRUE(subset)) {
            # Recalculate DESeq2 normalized counts and variance stabilizations
            # if the number of samples and/or genes change.
            if (isTRUE(recalculate)) {
                message("Recalculating DESeq2 normalizations.")
                dds <- .new.DESeqDataSet(se = rse)
                dds <- DESeq(dds)
                # Normalized counts.
                assays[["normalized"]] <- counts(dds, normalized = TRUE)
                # vst: variance-stabilizing transformation.
                if ("vst" %in% names(assays)) {
                    message("Applying variance-stabilizing transformation.")
                    assays[["vst"]] <-
                        assay(varianceStabilizingTransformation(dds))
                }
                # rlog: regularized log transformation.
                if ("rlog" %in% names(assays)) {
                    message("Applying rlog transformation.")
                    assays[["rlog"]] <- assay(rlog(dds))
                }
            } else {
                # Otherwise, remove previous calculations.
                message("Skipping DESeq2 normalizations.")
                assays[["normalized"]] <- NULL
                assays[["rlog"]] <- NULL
                assays[["vst"]] <- NULL
                assays[["fpkm"]] <- NULL
            }
        }

        # Row data -------------------------------------------------------------
        rowRanges <- rowRanges(rse)
        # Ensure factors get releveled, if necessary.
        # TODO Consider making this a function in basejump.
        if (nrow(rse) < nrow(x)) {
            mcols <- mcols(rowRanges)
            mcols <- DataFrame(lapply(
                X = mcols,
                FUN = function(x) {
                    if (is(x, "Rle")) {
                        x <- decode(x)
                        if (is.factor(x)) {
                            x <- droplevels(x)
                        }
                        Rle(x)
                    } else {
                        I(x)
                    }
                }
            ))
            mcols(rowRanges) <- mcols
        }

        # Column data ----------------------------------------------------------
        colData <- colData(rse)
        # Ensure factors get releveled, if necessary.
        if (ncol(rse) < ncol(x)) {
            colData <- colData %>%
                as.data.frame() %>%
                rownames_to_column() %>%
                mutate_if(is.character, as.factor) %>%
                mutate_if(is.factor, droplevels) %>%
                column_to_rownames() %>%
                as("DataFrame")
        }

        # Metadata -------------------------------------------------------------
        metadata <- metadata(rse)
        metadata[["subset"]] <- TRUE

        # Return ---------------------------------------------------------------
        .new.bcbioRNASeq(
            assays = assays,
            rowRanges = rowRanges,
            colData = colData,
            metadata = metadata
        )
    }
)
