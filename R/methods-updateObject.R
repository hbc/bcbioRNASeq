#' Update Object
#'
#' Update old objects created by the bcbioRNASeq package. The session
#' information metadata is preserved from the time when the bcbio data was
#' originally loaded into R.
#'
#' @section Upgrade bcbioRNADataSet to bcbioRNASeq object:
#' This method adds support for upgrading [bcbioRNADataSet] objects to the
#' latest [bcbioRNASeq] class version. This should be backwards compatible to
#' [bcbioRNASeq] version 0.0.26.
#'
#' @note The previous bcbioRnaseq package (note lowercase "c") must be
#'   reinstalled to load objects from versions <= 0.0.24.
#'
#' @rdname updateObject
#' @name updateObject
#' @author Michael Steinbaugh
#'
#' @importFrom BiocGenerics updateObject
#'
#' @inheritParams general
#'
#' @return [bcbioRNASeq].
#'
#' @examples
#' load(system.file(
#'     file.path("extdata", "bcb.rda"),
#'     package = "bcbioRNASeq"))
#'
#' metadata(bcb)[["version"]]
#'
#' updated <- updateObject(bcb)
#' metadata(updated)[["version"]]
#' metadata(updated)[["previousVersion"]]
NULL



# Constructors =================================================================
.updateObject.bcbioRNASeq <- function(object) {
    version <- metadata(object)[["version"]]
    assert_is_all_of(version, c("package_version", "numeric_version"))
    inform(paste("Upgrading from", version, "to", packageVersion))

    # Regenerate the bcbioRNASeq object
    se <- as(object, "SummarizedExperiment")
    to <- new("bcbioRNASeq", se)

    # Upgrade the bcbio slot
    bcbio <- bcbio(object)
    assert_is_all_of(bcbio, "SimpleList")
    slot(to, "bcbio") <- bcbio

    # Update the automatic metadata slots
    metadata(to)[["version"]] <- packageVersion
    metadata(to)[["previousVersion"]] <- metadata(object)[["version"]]
    metadata(to)[["upgradeDate"]] <- Sys.Date()

    # Version-specific modifications
    if (version <= package_version("0.0.26")) {
        # Remove GTF file, if present (too large)
        metadata(to)[["gtf"]] <- NULL
    }

    validObject(to)
    to
}



# Methods ======================================================================
#' @rdname updateObject
#' @export
setMethod(
    "updateObject",
    signature("bcbioRNASeq"),
    .updateObject.bcbioRNASeq)
