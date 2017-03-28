#' Import a bcbio project summary
#'
#' @author Michael Steinbaugh
#'
#' @import basejump
#' @import dplyr
#' @import readr
#'
#' @param bcbio bcbio run object
#'
#' @export
#' @examples
#' \dontrun{
#' import_summary(bcbio)
#' }
import_summary <- function(bcbio) {
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

    # Use the first entry in interesting groups to define the QC plot colors. If
    # empty, defaults to description.
    qc_color <- bcbio$intgroup[1]
    if (is.na(qc_color)) {
        qc_color = "description"
    }
    summary$qc_color <- metadata[rownames(summary), qc_color]

    return(summary)
}
