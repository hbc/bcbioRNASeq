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
#' loadRemoteData("http://bcbiornaseq.seq.cloud/f1000v1/bcb.rda", quiet = TRUE)
#' metadata(bcb)[["version"]]
#' updated <- updateObject(bcb)
#' metadata(updated)[["version"]]
#' metadata(updated)[["previousVersion"]]
NULL



# Constructors =================================================================
.updateObject.bcbioRNASeq <- function(object) {
    version <- metadata(object)[["version"]]
    assert_is_all_of(version, c("package_version", "numeric_version"))
    if (identical(version, packageVersion)) {
        inform("bcbioRNASeq object is already current")
        return(object)
    }

    inform(paste("Upgrading from", version, "to", packageVersion))

    # Regenerate the bcbioRNASeq object
    se <- as(object, "SummarizedExperiment")

    # Upgrade the metadata
    metadata <- metadata(se)
    assert_is_non_empty(names(metadata))

    # design
    if (is.null(metadata[["design"]])) {
        inform("Setting design as empty formula")
        metadata[["design"]] <- formula(~1)  # nolint
    }

    # ensemblVersion
    if (is.character(metadata[["ensemblVersion"]])) {
        inform("Setting ensemblVersion as NULL")
        metadata[["ensemblVersion"]] <- NULL
    } else if (
        !is.integer(metadata[["ensemblVersion"]]) &&
        is.numeric(metadata[["ensemblVersion"]])) {
        inform("Coercing ensemblVersion to integer")
        metadata[["ensemblVersion"]] <- as.integer(metadata[["ensemblVersion"]])
    }

    # gtf
    if ("gtf" %in% names(metadata)) {
        inform("Removing stashed GTF")
        metadata <- metadata[setdiff(names(metadata), "gtf")]
    }

    # lanes
    if (!is.integer(metadata[["lanes"]]) &&
        is.numeric(metadata[["lanes"]])) {
        inform("Setting lanes as integer")
        metadata[["lanes"]] <- as.integer(metadata[["lanes"]])
    }

    # programVersions
    if (!"programVersions" %in% names(metadata) &&
        "programs" %in% names(metadata)) {
        inform("Renaming programs to programVersions")
        metadata[["programVersions"]] <- metadata[["programs"]]
        metadata <- metadata[setdiff(names(metadata), "programs")]
    }

    # transformationLimit
    if (is.null(metadata[["transformationLimit"]])) {
        inform("Setting transformationLimit as Inf")
        metadata[["transformationLimit"]] <- Inf
    }

    # unannotatedGenes
    if (!"unannotatedGenes" %in% names(metadata) &&
        "missingGenes" %in% names(metadata)) {
        inform("Renaming missingGenes to unannotatedGenes")
        metadata[["unannotatedGenes"]] <- metadata[["missingGenes"]]
        metadata <- metadata[setdiff(names(metadata), "missingGenes")]
    }

    # Update the automatic metadata slots
    metadata[["version"]] <- packageVersion
    metadata[["previousVersion"]] <- metadata(object)[["version"]]
    metadata[["upgradeDate"]] <- Sys.Date()
    metadata(se) <- metadata

    # Upgrade the bcbio slot
    bcbio <- bcbio(object)
    assert_is_all_of(bcbio, "SimpleList")

    to <- new("bcbioRNASeq", se, bcbio = bcbio)
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
