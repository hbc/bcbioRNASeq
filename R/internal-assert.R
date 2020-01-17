## Check that the user input contains valid countsFromAbundance parameter used
## with tximport, for coercion to DESeqDataSet, DGEList objects.
## Updated 2020-01-17.
.assertHasValidCFA <- function(object) {
    assert(is(object, "SummarizedExperiment"))
    cfa <- metadata(object)[["countsFromAbundance"]]
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



## Updated 2020-01-17.
.isGeneLevel <- function(object) {
    identical(metadata(object)[["level"]], "genes")
}



## Updated 2020-01-17.
.isTranscriptLevel <- function(object) {
    identical(metadata(object)[["level"]], "transcripts")
}
