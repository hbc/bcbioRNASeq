# isGeneLevel ==================================================================
.isGeneLevel <- function(object) {
    identical(metadata(object)[["level"]], "genes")
}

.assertIsGeneLevel <- function(object) {
    validObject(object)
    assert_is_all_of(object, "bcbioRNASeq")
    assert_all_are_true(.isGeneLevel(object))
}



# isTranscriptLevel ============================================================
.isTranscriptLevel <- function(object) {
    identical(metadata(object)[["level"]], "transcripts")
}

.assertIsTranscriptLevel <- function(object) {
    validObject(object)
    assert_is_all_of(object, "bcbioRNASeq")
    assert_all_are_true(.isTranscriptLevel(object))
}
