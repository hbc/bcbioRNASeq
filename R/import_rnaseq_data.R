#' Import RNA-Seq data from a bcbio project
#'
#' @author Michael Steinbaugh
#'
#' @param bcbio bcbio run object
#'
#' @return List containing RNA-Seq data
#' @export
#'
#' @examples
#' \dontrun{
#' import_rnaseq_data(bcbio)
#' }
import_rnaseq_data <- function(bcbio) {
    # `.counts` = `featureCounts`
    annotated_combined.counts <-
        import_file(bcbio,
                    file = "annotated_combined.counts")
    combined.counts <-
        import_file(bcbio,
                    file = "combined.counts",
                    output = "matrix",
                    rownames = "id")

    # `.sf` = `sailfish`
    # Don't import the `combined.sf` file. Use the `tximport()` method on the
    # sailfish data in the sample folders instead.
    combined.gene.sf.tpm <-
        import_file(bcbio,
                    file = "combined.gene.sf.tpm",
                    output = "matrix",
                    rownames = "gene_id")
    combined.isoform.sf.tpm <-
        import_file(bcbio,
                    file = "combined.isoform.sf.tpm",
                    output = "matrix",
                    rownames = "id")

    return(list(
        annotated_combined.counts = annotated_combined.counts,
        combined.counts = combined.counts,
        combined.gene.sf.tpm = combined.gene.sf.tpm,
        combined.isoform.sf.tpm = combined.isoform.sf.tpm
    ))
}
