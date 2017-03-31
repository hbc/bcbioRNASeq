# Modified from bcbio-rnaseq qc-summary template
# https://goo.gl/cIv3ia

#' Import counts
#'
#' Import \code{bcbio-rnaseq} count data run using \code{tximport}
#'
#' @author Michael Steinbaugh
#' @author Rory Kirchner
#'
#' @import readr
#' @import stringr
#' @import tximport
#' @importFrom stats na.omit
#' @importFrom utils write.csv
#'
#' @param bcbio bcbio run object
#' @param type The type of software used to generate the abundances. Follows the
#'   conventions of \code{tximport()}.
#' @param pattern Apply pattern matching to sample names. This will override any
#'   parmeter set in `samples`.
#' @param samples Specify the names of samples in bcbio final directory to
#'   input. If \code{NULL} (default), all samples will be imported.
#'
#' @return txi \code{tximport} list object
#' @export
#'
#' @examples
#' \dontrun{
#' import_counts(bcbio, type = "salmon")
#' }
import_counts <- function(
    bcbio,
    type = "salmon",
    pattern = NULL,
    samples = NULL,
    save_tpm = FALSE) {
    check_bcbio_object(bcbio)
    if (!type %in% c("salmon", "sailfish")) {
        stop("unsupported counts input format")
    }

    # Sample name pattern matching. Run above `samples` to override.
    if (!is.null(pattern)) {
        samples <- stringr::str_subset(names(bcbio$sample_dirs),
                                       pattern = pattern)
        if (!length(samples)) {
            stop("pattern didn't match any samples")
        }
    }

    if (is.null(samples)) {
        samples <- names(bcbio$sample_dirs)
    } else {
        samples <- sort(unique(stats::na.omit(samples)))
    }

    # Draft support for salmon and sailfish file structure
    sample_files <- file.path(bcbio$sample_dirs[samples],
                              type, "quant", "quant.sf")
    names(sample_files) <- names(bcbio$sample_dirs[samples])
    if (!all(file.exists(sample_files))) {
        stop(paste(type, "count files do not exist"))
    }

    # Begin loading of selected counts
    message(paste("loading", type, "counts"))
    message(paste(names(sample_files), collapse = "\n"))
    # writeLines(sample_files)

    tx2gene <- import_file(bcbio, file = "tx2gene.csv", col_names = FALSE)

    # Import the counts
    # https://goo.gl/h6fm15
    # countsFromAbundance = c("no", "scaledTPM", "lengthScaledTPM")
    # Use `lengthScaledTPM` for `salmon` and `sailfish`
    txi <- tximport::tximport(files = sample_files,
                              type = type,
                              tx2gene = tx2gene,
                              reader = readr::read_tsv,
                              countsFromAbundance = "lengthScaledTPM")

    return(txi)
}
