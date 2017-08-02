#' Create `isomiRs` Object from `bcbio` Output
#'
#' Read `bcbio` sample information from YAML to get `isomiR` object.
#'
#' @rdname readSmallRNACounts
#' @author Lorena Patano
#' @keywords internal
#'
#' @param rna [bcbioRNADataSet].
.readSmallRNACounts <- function(rna) {
    meta <- metadata(rna)
    fns <- file.path(meta[["sampleDirs"]],
                     paste(names(meta[["sampleDirs"]]),
                           "mirbase-ready.counts",
                           sep = "-"))
    names(fns) <- names(meta[["sampleDirs"]])
    message("Reading miRNA count files")
    bcbio(rna, type = "isomirs") <- IsomirDataSeqFromFiles(
        files = fns[rownames(meta[["metadata"]])],
        coldata = meta[["metadata"]],
        design = ~meta[["interestingGroups"]])
    rna
}
