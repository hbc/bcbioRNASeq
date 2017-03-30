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
#' @import tximport
#' @importFrom utils write.csv
#'
#' @param bcbio bcbio run object
#' @param type The type of software used to generate the abundances. Follows the
#'   conventions of \code{tximport()}.
#' @param samples Specify the names of samples in bcbio final directory to
#'   input. If \code{NULL} (default), all samples will be imported.
#' @param save_tpm Whether to save transcripts per million data frame
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
    samples = NULL,
    save_tpm = FALSE) {
    if (!is.list(bcbio)) {
        stop("bcbio run object is required.")
    }
    if (!type %in% c("salmon", "sailfish")) {
        stop("Unsupported type.")
    }

    if (is.null(samples)) {
        samples <- names(bcbio$sample_dirs)
    }

    # Draft support for salmon and sailfish file structure
    sample_files <- file.path(bcbio$sample_dirs[samples],
                              type, "quant", "quant.sf")
    names(sample_files) <- names(bcbio$sample_dirs[samples])
    if (!all(file.exists(sample_files))) {
        stop(paste(type, "count files do not exist."))
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

    # Transcripts per million
    if (isTRUE(save_tpm)) {
        tpm <- txi$abundance
        assign("tpm", tpm, envir = parent.frame())
        save(tpm, file = "data/tpm.rda")
        # `write.csv()` required to write rownames
        utils::write.csv(tpm, file = "results/tpm.csv")
    }

    # Update bcbio object
    bcbio$sample_files <- sample_files
    # Assign sample_files into bcbio object
    assign("bcbio", bcbio, envir = parent.frame())
    save(bcbio, file = "data/bcbio.rda")

    save(txi, file = "data/txi.rda")
    return(txi)
}
