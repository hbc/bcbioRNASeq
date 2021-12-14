#' Assert that counts from abundance is valid
#'
#' @note Updated 2021-09-01.
#' @noRd
#'
#' @details
#' Check that the user input contains valid countsFromAbundance parameter used
#' with tximport, for coercion to DESeqDataSet, DGEList objects.
.assertHasValidCFA <- function(object) {
    assert(is(object, "SummarizedExperiment"))
    cfa <- metadata(object)[["countsFromAbundance"]]
    ## Note that featureCounts callers will return NULL here.
<<<<<<< HEAD
    if (is.null(cfa)) {
        return(TRUE)
    }
=======
    if (is.null(cfa)) return(TRUE)
>>>>>>> dc75f61043d7 (Improve error message handling)
    assert(
        isCharacter(cfa),
        isSubset(cfa, c("lengthScaledTPM", "no")),
        msg = sprintf(
            fmt = paste0(
                "Unsupported '%s' type: '%s'.\n",
                "Use either '%s' or '%s'. ",
                "See '%s' and '%s' documentation for details."
            ),
            "countsFromAbundance", cfa,
            "lengthScaledTPM", "no",
            "bcbioRNASeq()", "tximport()"
        )
    )
    TRUE
}



#' Assert that the bcbioRNASeq object was not generated using fast mode
#'
#' @note Updated 2021-09-01.
#' @noRd
.assertIsNotFastMode <- function(object) {
    assert(
        is(object, "bcbioRNASeq"),
        !.isFastMode(object),
        msg = sprintf(
            fmt = paste0(
                "This function does not support a %s object generated ",
                "with fast mode.\nRerun '%s' with '%s'."
            ),
            "bcbioRNASeq", "bcbioRNASeq()", "fast = FALSE"
        )
    )
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
        y = c("aligned", .deseqAssays)
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
    if (!hasLength(log)) {
        return(FALSE)
    }
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
