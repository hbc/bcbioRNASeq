# Default naming scheme is `template/final/YYYY-MM-DD_template`

#' Load bcbio run
#'
#' We recommend loading the \code{bcbio-nextgen} run as a remote connection over
#' \code{sshfs}. This requires setting \code{parent_dir} in the function.
#'
#' @author Michael Steinbaugh
#'
#' @param final_dir Path to final output directory. This path is set when
#'   running \code{bcbio_nextgen -w template}.
#' @param organism Organism (e.g. \code{hsapiens})
#' @param intgroup Character vector of interesting groups. First entry is used
#'   for plot colors during quality control (QC) analysis. Entire vector is used
#'   for PCA and heatmap QC functions.
#'
#' @return \code{bcbio-nextgen} run object
#' @export
load_bcbio_run <- function(
    final_dir = "final",
    organism,
    intgroup = "description") {
    if (!length(dir(final_dir))) {
        stop("final directory failed to load")
    }
    final_dir <- normalizePath(final_dir)
    message(paste(dir(final_dir), collapse = "\n"))

    # project_dir
    pattern <- ".*/(\\d{4}-\\d{2}-\\d{2})_([^/]+)$"
    match <- dir(final_dir, full.names = TRUE) %>%
        str_subset(pattern) %>%
        str_match(pattern)
    project_dir <- match[1]
    run_date <- match[2]
    template <- match[3]

    # Sample directories
    sample_dirs <- dir(final_dir, full.names = TRUE)
    names(sample_dirs) <- basename(sample_dirs)
    # Remove the nested `project_dir`
    sample_dirs <- sample_dirs[!grepl(basename(project_dir),
                                      names(sample_dirs))]

    # Data versions
    if (file.exists(file.path(project_dir, "data_versions.csv"))) {
        data_versions <- file.path(project_dir, "data_versions.csv") %>%
            read_csv(col_types = cols())
    } else {
        data_versions <- NULL
    }

    # Detect lane split FASTQ samples
    if (any(str_detect(dir(final_dir), "_L\\d+$"))) {
        lane_split <- TRUE
    } else {
        lane_split <- FALSE
    }

    run <- list(
        final_dir = final_dir,
        project_dir = project_dir,
        run_date = run_date,
        template = template,
        sample_dirs = sample_dirs,
        # User input
        organism = organism,  # Automate by parsing genome in YAML?
        intgroup = intgroup,
        # Fully automatic
        lane_split = lane_split,
        today_date = Sys.Date(),
        session_info = sessionInfo(),
        data_versions = data_versions,
        hpc = detect_hpc(),
        wd = getwd()
    )

    check_run(run)
    return(run)
}
