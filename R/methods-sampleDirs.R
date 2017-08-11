#' Sample Directories
#'
#' This method will be used to access folders where sample information is kept.
#'
#' @rdname sampleDirs
#' @name sampleDirs
#'
#' @return Named character vector containing sample directory paths.
#'
#' @examples
#' data(bcb)
#' sampleDirs(bcb) %>% basename
NULL



# Methods ====
#' @rdname sampleDirs
#' @export
setMethod("sampleDirs", "bcbioRNADataSet", function(object) {
    metadata(object)[["sampleDirs"]]
})
