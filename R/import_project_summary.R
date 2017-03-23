#' Import a bcbio project summary
#'
#' @author Michael Steinbaugh
#'
#' @import basejump
#' @import dplyr
#' @import readr
#'
#' @param bcbio bcbio run object
#' @param metadata bcbio project metadata
#'
#' @export
#' @examples
#' \dontrun{
#' import_project_summary(bcbio)
#' }
import_project_summary <- function(bcbio, metadata) {
    if (!is.data.frame(metadata)) {
        stop("A metadata data frame is required.")
    }
    df <- file.path(bcbio$project_dir,
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

    # Set the group, used for plots
    df$group <- metadata[rownames(df), "group"]

    return(df)
}
