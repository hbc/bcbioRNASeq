.assertIsGeneLevel <- function(object) {
    validObject(object)
    assert_is_all_of(object, "bcbioRNASeq")
    assert_all_are_true(.isGeneLevel(object))
}

.isGeneLevel <- function(object) {
    identical(metadata(object)[["level"]], "genes")
}



.assertIsTranscriptLevel <- function(object) {
    validObject(object)
    assert_is_all_of(object, "bcbioRNASeq")
    assert_all_are_true(.isTranscriptLevel(object))
}

.isTranscriptLevel <- function(object) {
    identical(metadata(object)[["level"]], "transcripts")
}
