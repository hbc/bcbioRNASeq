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
#' @param summary Summary data frame
#' @param ... Passthrough to \code{tximport()}
#'
#' @return txi \code{tximport} list object
#' @export
#'
#' @examples
#' \dontrun{
#' import_sailfish(bcbio,
#'                 type = "salmon",
#'                 summary = summary)
#' }
import_counts <- function(
    bcbio,
    type = "salmon",
    summary,
    ...) {
    if (!is.list(bcbio)) {
        stop("bcbio run object is required.")
    }
    if (!type %in% c("salmon", "sailfish")) {
        stop("Unsupported type.")
    }
    if (!is.data.frame(summary)) {
        stop("bcbio run summary data frame is required.")
    }

    # tx2gene
    if (!file.exists(file.path(bcbio$project_dir, "tx2gene.csv"))) {
        stop("tx2gene.csv file not found.")
    }
    tx2gene <- file.path(bcbio$project_dir, "tx2gene.csv") %>%
        utils::read.csv(header = FALSE)

    # Parse the HPC run
    files <- dir(bcbio$final_dir) %>%
        .[. %in% summary$description] %>%
        file.path(bcbio$final_dir,
                  .,
                  type,
                  "quant",
                  "quant.sf") %>%
        sort
    if (!length(files)) {
        stop(paste("No", type, "files were found."))
    }

    names(files) <- sort(summary$description)

    # Import the counts
    # https://goo.gl/h6fm15
    # countsFromAbundance = c("no", "scaledTPM", "lengthScaledTPM")
    txi <- tximport::tximport(files,
                              type = type,
                              tx2gene = tx2gene,
                              reader = readr::read_tsv,
                              countsFromAbundance = "lengthScaledTPM",
                              ...)

    save(txi, file = "data/txi.rda")
    return(txi)
}



#' Import sailfish counts
#'
#' @keywords internal
#' @param ... Passthrough to \code{import_counts()}
#' @export
import_sailfish <- function(...) {
    import_counts(..., type = "sailfish")
}

#' Import salmon counts
#'
#' @keywords internal
#' @param ... Passthrough to \code{import_counts()}
#' @export
import_salmon <- function(...) {
    import_counts(..., type = "salmon")
}
