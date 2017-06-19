#' Create isomiRs object from bcbio output
#'
#' Read bcbio sample information from YAML to get isomiR object.
#'
#' @rdname read_smallrna_counts
#' @keywords internal
#'
#' @author Lorena Patano
#'
#' @param rna [bcbioRnaDataSet].
.read_smallrna_counts <- function(rna) {
    # [TODO] Better way to handle sample_dirs than by piping in via metadata?
    run <- metadata(rna)
    fns <- file.path(run$sample_dirs,
                     paste(names(run$sample_dirs),
                           "mirbase-ready.counts",
                           sep = "-"))
    names(fns) <- names(run$sample_dirs)
    message("Reading miRNA count files...")
    bcbio(rna, type = "isomirs") <- IsomirDataSeqFromFiles(
        files = fns[rownames(run$metadata)],
        coldata = run$metadata,
        design = ~run$groups_of_interest)
    rna
}
