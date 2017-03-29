#' Import a bcbio project summary
#'
#' @author Michael Steinbaugh
#'
#' @import basejump
#' @import dplyr
#' @import readr
#' @importFrom utils write.csv
#'
#' @param bcbio bcbio run object
#' @param save Save data frame
#'
#' @export
#' @examples
#' \dontrun{
#' import_summary(bcbio)
#' }
import_summary <- function(bcbio, save = FALSE) {
    summary <- file.path(bcbio$project_dir,
                         "project-summary.csv") %>%
        readr::read_csv(., col_types = readr::cols()) %>%
        basejump::setNamesSnake(.) %>%
        # Remove NA only columns
        .[, colSums(!is.na(.)) > 0] %>%
        # Sort by description
        dplyr::select_(.dots = c("description",
                                 setdiff(sort(names(.)),
                                         "description"))) %>%
        dplyr::arrange_(.dots = "description") %>%
        set_rownames("description")

    metadata <- import_metadata(bcbio)

    if (!identical(rownames(summary), rownames(metadata))) {
        stop("Summary and metadata rownames don't match.")
    }

    # Use the first entry in interesting groups to define the QC plot colors. If
    # empty, defaults to description.
    color <- bcbio$intgroup[1]
    if (is.na(color)) {
        color = "description"
    }
    summary$qc_color <- metadata[rownames(summary), color]

    if (isTRUE(save)) {
        save(summary, file = "data/summary.rda")
        utils::write.csv(summary, file = "results/summary.csv")
    }

    return(summary)
}
