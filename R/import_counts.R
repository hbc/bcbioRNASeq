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
#' @importFrom utils read.csv
#'
#' @param bcbio bcbio run object
#' @param type The type of software used to generate the abundances. Follows the
#'   conventions of \code{tximport()}.
#' @param ... Passthrough to \code{tximport()}
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
    ...) {
    if (!is.list(bcbio)) {
        stop("bcbio run object is required.")
    }
    if (!type %in% c("salmon", "sailfish")) {
        stop("Unsupported type.")
    }

    summary <- import_summary(bcbio)

    # tx2gene
    if (!file.exists(file.path(bcbio$project_dir, "tx2gene.csv"))) {
        stop("tx2gene.csv file not found.")
    }
    tx2gene <- file.path(bcbio$project_dir, "tx2gene.csv") %>%
        utils::read.csv(header = FALSE)

    # Parse the bcbio RNA-Seq run
    # bcbio names the final sample folders by `description`
    files <- dir(bcbio$final_dir) %>%
        .[. %in% summary$description] %>%
        file.path(bcbio$final_dir,
                  .,
                  type,
                  "quant",
                  "quant.sf") %>%
        sort
    if (!length(files)) {
        stop(paste(type, "files failed to load."))
    }

    names(files) <- sort(summary$description)

    # Import the counts
    # https://goo.gl/h6fm15
    # countsFromAbundance = c("no", "scaledTPM", "lengthScaledTPM")
    # Use `lengthScaledTPM` for `salmon` and `sailfish`
    txi <- tximport::tximport(files,
                              type = type,
                              tx2gene = tx2gene,
                              reader = readr::read_tsv,
                              countsFromAbundance = "lengthScaledTPM",
                              ...)

    # Transcripts per million
    tpm <- txi$abundance
    assign("tpm", tpm, envir = parent.frame())

    if (dir.exists("data")) {
        save(tpm, file = "data/tpm.rda")
        save(txi, file = "data/txi.rda")
    }

    if (dir.exists("results")) {
        write.csv(tpm, file = "results/tpm.csv")
    }

    return(txi)
}
