#' Update Object
#'
#' Compatible with old versions created by bcbioRNASeq package.
#'
#' @rdname updateObject
#' @name updateObject
#' @author Michael Steinbaugh
#'
#' @importFrom BiocGenerics updateObject
#'
#' @note The previous bcbioRnaseq package (note lowercase "c") must be
#'   reinstalled to load objects from versions <= 0.0.24.
#'
#' @return [bcbioRNASeq].
#'
#' @examples
#' # FIXME Needs working example
NULL



# Constructors =================================================================
#' @importFrom S4Vectors metadata
.updateObject.bcbioRNADataSet <- function(object) {
    version <- metadata(object)[["version"]]
    assert_is_all_of(version, c("package_version", "numeric_version"))

    inform(paste(
        paste("Upgrading from", version, "to", packageVersion),
        paste("Existing metadata:", toString(names(metadata(object)))),
        sep = "\n"
    ))

    # Regenerate the bcbioRNASeq object
    se <- as(object, "SummarizedExperiment")
    to <- new("bcbioRNASeq", se)
    bcbio(to) <- bcbio(object)
    validObject(to)

    # Update the automatic metadata slots
    metadata(to)[["version"]] <- packageVersion
    metadata(to)[["originalVersion"]] <- metadata(object)[["version"]]
    metadata(to)[["upgradeDate"]] <- Sys.Date()

    # Version-specific modifications
    if (version <= package_version("0.0.26")) {
        # Remove GTF file, if present (too large)
        metadata(to)[["gtf"]] <- NULL
    }

    to
}



# Methods ======================================================================
#' @rdname updateObject
#'
#' @section Upgrade bcbioRNADataSet to bcbioRNASeq object:
#' This method adds support for upgrading [bcbioRNADataSet] objects to the
#' latest [bcbioRNASeq] class version. This should be backwards compatible to
#' [bcbioRNASeq] version 0.0.26. Previous objects saved using `bcbioRnaseq`
#' (note case) will likely fail to load with newer versions of the package.
#'
#' @export
setMethod(
    "updateObject",
    signature("bcbioRNADataSet"),
    .updateObject.bcbioRNADataSet)
