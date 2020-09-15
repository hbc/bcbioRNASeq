#' Assert that counts from abundance is valid
#'
#' @note Updated 2020-09-15.
#' @noRd
#'
#' @details
#' Check that the user input contains valid countsFromAbundance parameter used
#' with tximport, for coercion to DESeqDataSet, DGEList objects.
.assertHasValidCFA <- function(object) {
    assert(is(object, "SummarizedExperiment"))
    cfa <- metadata(object)[["countsFromAbundance"]]
    ## Note that featureCounts callers will return NULL here.
    if (is.null(cfa)) return(TRUE)
    assert(isCharacter(cfa))
    if (!isSubset(cfa, c("lengthScaledTPM", "no"))) {
        stop(
            "Unsupported 'countsFromAbundance' type: ", cfa, ".\n",
            "Use either 'lengthScaledTPM' or 'no'. ",
            "See `bcbioRNASeq()` and `tximport()` documentation for details."
        )
    }
    TRUE
}



#' Assert that the bcbioRNASeq object was not generated using fast mode
#'
#' @note Updated 2020-09-15.
#' @noRd
.assertIsNotFastMode <- function(object) {
    assert(is(object, "bcbioRNASeq"))
    if (.isFastMode(object)) {
        stop(paste0(
            "This function does not support a bcbioRNASeq object generated ",
            "with fast mode. Rerun 'bcbioRNASeq()' with 'fast = FALSE'."
        ))
    }
    TRUE
}



#' Was the bcbioRNASeq object generated using fast mode?
#'
#' @note Updated 2020-09-15.
#' @noRd
.isFastMode <- function(object) {
    assert(is(object, "bcbioRNASeq"))
    areDisjointSets(
        x = assayNames(object),
        y = c("aligned", "fpkm", "normalized", "vst")
    )
}



#' Is bcbio-nextgen output the fastrnaseq pipeline?
#'
#' @note Updated 2020-09-15.
#' @noRd
#'
#' @details
#' Check the `bcbio-nextgen.log` file to see if fastrnaseq pipeline was run.
.isFastPipeline <- function(log) {
    if (!hasLength(log)) return(FALSE)
    any(grepl(pattern = "fastrnaseq", x = log, fixed = TRUE))
}



#' Is the bcbioRNASeq object at gene level?
#'
#' @note Updated 2020-01-17.
#' @noRd
.isGeneLevel <- function(object) {
    identical(metadata(object)[["level"]], "genes")
}



#' Is the bcbioRNASeq object at transcript level?
#'
#' @note Updated 2020-01-17.
#' @noRd
.isTranscriptLevel <- function(object) {
    identical(metadata(object)[["level"]], "transcripts")
}
