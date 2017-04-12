#' @rdname run_import
#' @description Import project quality control summary statistics
#' @return Summary data frame
#' @export
import_qc_summary <- function(run, save = FALSE) {
    check_run(run)
    summary <- file.path(run$project_dir,
                         "project-summary.csv") %>%
        read_csv(col_types = cols()) %>%
        set_names_snake %>%
        # Remove NA only columns
        .[, colSums(!is.na(.)) > 0] %>%
        # Sort by description
        select_(.dots = c("description",
                          setdiff(sort(names(.)),
                                  "description"))) %>%
        arrange_(.dots = "description")

    metadata <- import_metadata(run)

    if (!identical(summary$description, metadata$description)) {
        stop("summary and metadata descriptions don't match")
    }

    # Use the first entry in interesting groups to define the QC colors. If
    # empty, defaults to description.
    color <- run$intgroup[1]
    if (is.na(color)) {
        color <- "description"
    }
    summary$qc_color <- metadata[[color]]

    if (isTRUE(save)) {
        save(summary, file = "data/qc_summary.rda")
        write_csv(summary, "results/qc_summary.csv")
    }

    return(summary)
}
