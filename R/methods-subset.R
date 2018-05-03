#' Bracket-Based Subsetting
#'
#' Extract genes by row and samples by column from a `bcbioRNASeq` object. The
#' internal `DESeqDataSet` and count transformations are rescaled automatically.
#' DESeq2 transformations can be disabled on large subset operations by setting
#' `transform = FALSE`.
#'
#' @name subset
#' @family S4 Class Definition
#' @author Lorena Pantano, Michael Steinbaugh
#'
#' @inheritParams base::`[`
#' @inheritParams general
#'
#' @return `bcbioRNASeq`.
#'
#' @seealso `help("[", "base")`.
#'
#' @examples
#' # Minimum of 100 genes, 2 samples
#' genes <- head(rownames(bcb_small), 100L)
#' head(genes)
#' samples <- head(colnames(bcb_small), 2L)
#' head(samples)
#'
#' # Subset by sample name
#' bcb_small[, samples]
#'
#' # Subset by gene list
#' bcb_small[genes, ]
#'
#' # Subset by both genes and samples
#' subset <- bcb_small[genes, samples]
#' print(subset)
#' assayNames(subset)
NULL



# Methods ======================================================================
#' @rdname subset
#' @export
setMethod(
    "[",
    signature(
        x = "bcbioRNASeq",
        i = "ANY",
        j = "ANY",
        drop = "ANY"
    ),
    function(x, i, j, drop = FALSE) {
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

        # Early return if dimensions are unmodified
        if (identical(dim(x), c(length(i), length(j)))) {
            return(x)
        }

        # Regenerate RangedSummarizedExperiment
        rse <- as(x, "RangedSummarizedExperiment")
        rse <- rse[i, j, drop = drop]

        # Assays ===============================================================
        assays <- assays(rse)

        # Update DESeq2 transformations, if they are defined
        if (any(c("rlog", "vst") %in% names(assays))) {
            message("Updating variance stabilizations")
            dds <- .regenerateDESeqDataSet(rse)
            message("Applying rlog transformation")
            assays[["rlog"]] <- assay(rlog(dds))
            message("Applying variance stabilizing transformation")
            assays[["vst"]] <- assay(varianceStabilizingTransformation(dds))
        }

        # Column data ==========================================================
        # Ensure factors get releveled
        colData <- colData(rse) %>%
            as.data.frame() %>%
            rownames_to_column() %>%
            mutate_if(is.character, as.factor) %>%
            mutate_if(is.factor, droplevels) %>%
            column_to_rownames() %>%
            as("DataFrame")

        # Metadata =============================================================
        metadata <- metadata(rse)
        metadata[["subset"]] <- TRUE
        # Update version, if necessary
        if (!identical(metadata[["version"]], packageVersion)) {
            metadata[["originalVersion"]] <- metadata[["version"]]
            metadata[["version"]] <- packageVersion
        }

        # Return ===============================================================
        .new.bcbioRNASeq(
            assays = assays,
            rowRanges = rowRanges(rse),
            colData = colData,
            metadata = metadata
        )
    }
)
