.assertIsGeneLevel <- function(object) {
    validObject(object)
    assert_is_all_of(object, "bcbioRNASeq")
    assert_all_are_true(.isGeneLevel)
}

.assertIsTranscriptLevel <- function(object) {
    validObject(object)
    assert_is_all_of(object, "bcbioRNASeq")
    assert_all_are_true(.isTranscriptLevel)
}

.isGeneLevel <- function(object) {
    identical(metadata(object)[["level"]], "genes")
}

.isTranscriptLevel <- function(object) {
    identical(metadata(object)[["level"]], "transcripts")
}
