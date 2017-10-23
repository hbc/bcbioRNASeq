#' Sample Directories
#'
#' This method will be used to access folders where sample information is kept.
#'
#' @rdname sampleDirs
#' @name sampleDirs
#'
#' @importFrom basejump sampleDirs
#'
#' @inheritParams AllGenerics
#'
#' @return Named character vector containing sample directory paths.
#'
#' @examples
#' data(bcb)
#' sampleDirs(bcb) %>% basename()
NULL



# Methods ====
#' @rdname sampleDirs
#' @importFrom S4Vectors metadata
#' @export
setMethod(
    "sampleDirs",
    signature("bcbioRNASeq"),
    function(object) {
        metadata(object)[["sampleDirs"]]
    })
