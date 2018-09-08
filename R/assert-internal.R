.assertIsGeneLevel <- function(object) {
    validObject(object)
    assert_is_all_of(object, "bcbioRNASeq")
    assert_are_identical(
        x = metadata(object)[["level"]],
        y = "genes"
    )
}

.assertIsTranscriptLevel <- function(object) {
    validObject(object)
    assert_is_all_of(object, "bcbioRNASeq")
    assert_are_identical(
        x = metadata(object)[["level"]],
        y = "transcripts"
    )
}
