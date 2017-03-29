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
    save_tpm = TRUE) {
    if (!is.list(bcbio)) {
        stop("bcbio run object is required.")
    }
    if (!type %in% c("salmon", "sailfish")) {
        stop("Unsupported type.")
    }

    summary <- import_summary(bcbio)
    tx2gene <- import_file(bcbio, file = "tx2gene.csv", col_names = FALSE)

    # Parse the bcbio RNA-Seq run
    # bcbio names the final sample folders by `description`
    sample_dirs <- dir(bcbio$final_dir) %>%
        .[. %in% summary$description] %>% sort %>% unique
    if (!length(sample_dirs)) {
        stop("No sample directories matched the summary report.")
    }

    sample_files <- file.path(bcbio$final_dir, sample_dirs,
                              type, "quant", "quant.sf")
    names(sample_files) <- sample_dirs
    if (!all(file.exists(sample_files))) {
        stop(paste(type, "count files do not exist."))
    }

    message(paste("Loading", type, "counts..."))
    message(paste(sample_dirs, collapse = "\n"))
    # writeLines(sample_dirs)

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
        if (dir.exists("data")) {
            save(tpm, file = "data/tpm.rda")
        }
        if (dir.exists("results")) {
            # `write.csv()` required to write rownames
            utils::write.csv(tpm, file = "results/tpm.csv")
        }
    }

    if (dir.exists("data")) {
        save(txi, file = "data/txi.rda")
    }

    return(txi)
}
