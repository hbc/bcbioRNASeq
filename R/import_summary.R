#' @rdname import
#' @description Import project summary statistics
#'
#' @import dplyr
#' @import readr
#' @importFrom utils write.csv
#'
#' @param bcbio bcbio run object
#' @param save Save data frame
#'
#' @return Summary data frame
#' @export
import_summary <- function(bcbio, save = FALSE) {
    check_bcbio(bcbio)
    summary <- file.path(bcbio$project_dir,
                         "project-summary.csv") %>%
        read_csv(col_types = cols()) %>%
        set_names_snake %>%
        # Remove NA only columns
        .[, colSums(!is.na(.)) > 0] %>%
        # Sort by description
        select_(.dots = c("description",
                          setdiff(sort(names(.)),
                                  "description"))) %>%
        arrange_(.dots = "description") %>%
        set_rownames("description")

    metadata <- import_metadata(bcbio)

    if (!identical(rownames(summary), rownames(metadata))) {
        stop("summary and metadata rownames don't match")
    }

    # Use the first entry in interesting groups to define the QC colors. If
    # empty, defaults to description.
    color <- bcbio$intgroup[1]
    if (is.na(color)) {
        color = "description"
    }
    summary$qc_color <- metadata[rownames(summary), color]

    if (isTRUE(save)) {
        save(summary, file = "data/summary.rda")
        write.csv(summary, file = "results/summary.csv")
    }

    return(summary)
}
