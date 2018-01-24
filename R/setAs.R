#' Coerce Object
#'
#' @rdname coerce
#' @name coerce
#' @author Michael Steinbaugh
#'
#' @param from Class for which the coerce method will perform coercion.
#'
#' @seealso `help(topic = "coerce", package = "methods")`.
#'
#' @examples
#' \dontrun{
#' # Legacy bcbioRNADataSet class
#' new <- as(old, "bcbioRNASeq")
#' }
NULL



# Constructors =================================================================
#' Coerce Legacy bcbio Object to bcbioRNASeq class
#'
#' Compatible with old versions created by bcbioRNASeq package.
#' The previous bcbioRnaseq package (note lowercase "c") must be reinstalled
#' to load objects from versions <= 0.0.24.
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @importFrom S4Vectors metadata
#'
#' @return [bcbioRNASeq].
.coerceLegacy <- function(from) {
    # Check for version
    version <- metadata(from)[["version"]]
    if (is.null(version)) {
        abort(paste(
            "Unknown bcbio object version.",
            "Please reload with 'loadRNASeq()'."
        ))
    }
    inform(paste(
        paste("Upgrading from", version, "to", packageVersion),
        paste("Existing metadata:", toString(names(metadata(from)))),
        sep = "\n"
    ))

    # Regenerate the bcbioRNASeq object
    se <- as(from, "SummarizedExperiment")
    to <- new("bcbioRNASeq", se)
    bcbio(to) <- bcbio(from)
    validObject(to)

    # Update the automatic metadata slots
    metadata(to)[["version"]] <- packageVersion
    metadata(to)[["originalVersion"]] <- metadata(from)[["version"]]
    metadata(to)[["upgradeDate"]] <- Sys.Date()

    # Version-specific modifications
    if (version <= package_version("0.0.26")) {
        # Remove GTF file, if present (too large)
        metadata(to)[["gtf"]] <- NULL
    }

    to
}



.coerceToSummarizedExperiment <- function(from) {
    to <- new("SummarizedExperiment")
    slot(to, "colData") <- slot(from, "colData")
    slot(to, "assays") <- slot(from, "assays")
    slot(to, "NAMES") <- slot(from, "NAMES")
    slot(to, "elementMetadata") <- slot(from, "elementMetadata")
    slot(to, "metadata") <- slot(from, "metadata")
    validObject(to)
    to
}



# Methods ======================================================================
#' @rdname coerce
#' @name upgrade-bcbioRNASeq
#' @section Upgrade bcbioRNASeq to current version:
#' This method adds support for upgrading [bcbioRNADataSet] objects to the
#' latest [bcbioRNASeq] class version. This should be backwards compatible to
#' [bcbioRNASeq] version 0.0.26. Previous objects saved using `bcbioRnaseq`
#' (note case) will likely fail to load with newer versions of the package.
setAs(
    from = "bcbioRNADataSet",
    to = "bcbioRNASeq",
    .coerceLegacy)



#' @rdname coerce
#' @name coerce-bcbioRNASeq-SummarizedExperiment
#' @section bcbioRNASeq to SummarizedExperiment:
#' Since [bcbioRNASeq] is an extension of [SummarizedExperiment], this
#' coercion method is very simple. Here we're simply dropping our `@bcbio` slot,
#' which contains raw cellular barcodes and other bcbio-specific metadata.
setAs(
    from = "bcbioRNASeq",
    to = "SummarizedExperiment",
    .coerceToSummarizedExperiment)
