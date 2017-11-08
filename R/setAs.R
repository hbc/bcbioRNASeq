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
#' # bcbioRNADataSet
#' newClass <- as(oldClass, "bcbioRNASeq")
#' }
NULL



# Constructors ====
#' Coerce Legacy bcbio Object to `bcbioRNASeq` class
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
#' @importFrom utils packageVersion
#'
#' @return [bcbioRNASeq].
.coerceLegacy <- function(from) {
    # Check for version
    version <- metadata(from)[["version"]]
    if (is.null(version)) {
        stop(paste(
            "Unknown bcbio object version.",
            "Please reload with 'loadRNASeq()'."
        ), call. = FALSE)
    }
    message(paste(
        "Upgrading from",
        version,
        "to",
        packageVersion("bcbioRNASeq")
    ))
    message(paste("Existing metadata:", toString(names(metadata(from)))))

    assays <- assays(from)

    rowData <- rowData(from)
    rownames(rowData) <- slot(from, "NAMES")

    colData <- colData(from)

    metadata <- metadata(from)
    metadata[["originalVersion"]] <- metadata[["version"]]
    metadata[["version"]] <- packageVersion("bcbioRNASeq")
    metadata[["upgradeDate"]] <- Sys.Date()

    # Version-specific modifications ====
    if (version <= package_version("0.0.26")) {
        bcbio <- slot(from, "callers")
        # Remove GTF file, if present (too large)
        metadata[["gtf"]] <- NULL
    } else {
        bcbio <- slot(from, "bcbio")
    }

    se <- SummarizedExperiment(
        assays = assays,
        rowData = rowData,
        colData = colData,
        metadata = metadata)

    # Return updated object ====
    to <- new("bcbioRNASeq", se)
    slot(to, "bcbio") <- bcbio
    validObject(to)
    to
}



# Methods ====
#' @rdname coerce
#' @name upgrade-bcbioRNASeq
#' @section Upgrade `bcbioRNASeq` to current version:
#' This method adds support for upgrading [bcbioRNADataSet] objects to the
#' latest [bcbioRNASeq] class version. This should be backwards compatible to
#' [bcbioRNASeq] version 0.0.26. Previous objects saved using `bcbioRnaseq`
#' (note case) will likely fail to load with newer versions of the package.
setAs(
    "bcbioRNADataSet",
    signature("bcbioRNASeq"),
    .coerceLegacy)
