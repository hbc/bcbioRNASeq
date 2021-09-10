#' Extract or replace parts of an object
#'
#' Extract genes by row and samples by column.
#'
#' Internal count transformations are rescaled automatically, if defined. DESeq2
#' transformations will only be updated when `recalculate = TRUE`.
#'
#' @name extract
#' @author Michael Steinbaugh, Lorena Pantano
#' @inherit base::Extract params references
#' @note Updated 2020-12-21.
#'
#' @inheritParams AcidRoxygen::params
#' @param recalculate `logical(1)`.
#'   Recalculate DESeq2 normalized counts and variance-stabilizing
#'   transformations defined in `assays`.\cr
#'   Recommended by default, but can take a long time for large datasets.\cr\
#'   If `FALSE`, these assays will be removed automatically: `normalized`,
#'   `rlog`, `vst`.
#'
#' @return `bcbioRNASeq`.
#'
#' @examples
#' ## bcbioRNASeq ====
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
#' ## Note that `normalized`, `rlog`, and `vst` assays will be removed.
#' x <- object[, samples, recalculate = FALSE]
#' print(x)
#' names(assays(x))
NULL



## Updated 2020-12-21.
`extract,bcbioRNASeq` <-  # nolint
    function(
        x, i, j,
        drop = FALSE,
        recalculate = TRUE
    ) {
        validObject(x)
        assert(
            identical(drop, FALSE),
            isFlag(recalculate)
        )
        ## Genes (rows).
        if (missing(i)) {
            i <- seq_len(nrow(x))
        }
        ## Require at least 50 genes.
        assert(isInRange(length(i), lower = 50L, upper = Inf))
        ## Samples (columns).
        if (missing(j)) {
            j <- seq_len(ncol(x))
        }
        ## Require at least 2 samples.
        assert(isInRange(length(j), lower = 2L, upper = Inf))
        ## Determine whether we should stash subset in metadata.
        if (identical(x = dim(x), y = c(length(i), length(j)))) {
            subset <- FALSE
        } else {
            subset <- TRUE
        }
        ## Regenerate RangedSummarizedExperiment.
        rse <- as(x, "RangedSummarizedExperiment")
        rse <- rse[i, j, drop = FALSE]
        ## Early return original object, if unmodified.
        if (identical(assay(rse), assay(x))) {
            return(x)
        }
        ## Assays --------------------------------------------------------------
        assays <- assays(rse)
        if (isTRUE(subset)) {
            ## Recalculate DESeq2 normalized counts and variance stabilizations
            ## if the number of samples and/or genes change.
            if (isTRUE(recalculate)) {
                alert(sprintf(
                    "Recalculating {.pkg %s} normalizations.", "DESeq2"
                ))
                dds <- `new,DESeqDataSet`(se = rse)
                dds <- DESeq(dds)
                ## Normalized counts.
                assays[["normalized"]] <- counts(dds, normalized = TRUE)
                ## vst: variance-stabilizing transformation.
                if ("vst" %in% names(assays)) {
                    alert("{.fun varianceStabilizingTransformation}")
                    assays[["vst"]] <-
                        assay(varianceStabilizingTransformation(dds))
                }
                ## rlog: regularized log transformation.
                ## v0.3.22: We're no longer allowing automatic rlog calculation
                ## during the `bcbioRNASeq()` call, because it's often too slow.
                ## However, we're keeping support here for legacy objects.
                if ("rlog" %in% names(assays)) {
                    alert(sprintf("{.fun %s}", "rlog"))
                    assays[["rlog"]] <- assay(rlog(dds))
                }
            } else {
                ## Otherwise, remove previous calculations.
                alertWarning(sprintf(
                    "Skipping {.fun %s} normalizations.", "DESeq2"
                ))
                assays[["normalized"]] <- NULL
                assays[["rlog"]] <- NULL
                assays[["vst"]] <- NULL
                assays[["fpkm"]] <- NULL
            }
        }
        ## Drop any `NULL` items in assays.
        assays <- Filter(f = Negate(is.null), x = assays)
        assays(rse) <- assays
        ## Metadata ------------------------------------------------------------
        metadata <- metadata(rse)
        if (isTRUE(subset)) {
            metadata[["subset"]] <- TRUE
        }
        ## Drop any `NULL` items in metadata.
        metadata <- Filter(f = Negate(is.null), x = metadata)
        metadata(rse) <- metadata
        ## Return --------------------------------------------------------------
        rse <- droplevels(rse)
        bcb <- new(Class = "bcbioRNASeq", rse)
        validObject(bcb)
        bcb
    }



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
    definition = `extract,bcbioRNASeq`
)
