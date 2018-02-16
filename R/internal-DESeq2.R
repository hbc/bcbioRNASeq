#' `DESeqResults` Contrast Name
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @param object [DESeqResults].
#'
#' @return Contrast name string.
#'
#' @examples
#' load(system.file(
#'     file.path("extdata", "res.rda"),
#'     package = "bcbioRNASeq"))
#' .contrastName.DESeqResults(res)
.contrastName.DESeqResults <- function(object) {  # nolint
    assert_is_all_of(object, "DESeqResults")
    contrast <- mcols(object)[2L, 2L]
    assert_is_character(contrast)
    contrast %>%
        gsub("^.*:\\s", "", .) %>%
        gsub("_", " ", .) %>%
        # Improve appearance for difference of differences
        gsub("\\+", " \\+\n    ", .)
}
