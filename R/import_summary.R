#' Import a bcbio project summary
#'
#' @author Michael Steinbaugh
#'
#' @import basejump
#' @import dplyr
#' @import readr
#'
#' @param bcbio bcbio run object
#' @param metadata Metadata data frame
#'
#' @export
#' @examples
#' \dontrun{
#' import_summary(bcbio)
#' }
import_summary <- function(bcbio, metadata) {
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

    summary$intgroup <- metadata[rownames(summary), "intgroup"]
    return(summary)
}
