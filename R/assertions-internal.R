.isGeneLevel <- function(object) {
    identical(metadata(object)[["level"]], "genes")
}



.isTranscriptLevel <- function(object) {
    identical(metadata(object)[["level"]], "transcripts")
}
