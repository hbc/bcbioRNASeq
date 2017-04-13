#' Load bcbio run
#'
#' We recommend loading the \code{bcbio-nextgen} run as a remote connection over
#' \code{sshfs}. This requires setting \code{parent_dir} in the function.
#'
#' @author Michael Steinbaugh
#'
#' @param template Template name (YAML filename, without the extension)
#'   originally called by \code{bcbio_nextgen.py}
#' @param organism Organism (e.g. \code{hsapiens})
#' @param intgroup Interesting groups. Defaults to \code{description} if not
#'   set. First entry is used for plot colors during quality control analysis.
#'   Entire vector is used for PCA and heatmap QC functions.
#' @param parent_dir Parent directory containing the run folder
#' @param config_dir Set the config output directory, if non-standard
#' @param final_dir Set the final output directory, if non-standard
#'
#' @return \code{bcbio-nextgen} run object
#' @export
load_bcbio_run <- function(
    template,
    organism,
    intgroup = NULL,
    parent_dir = getwd(),
    config_dir = NULL,
    final_dir = NULL) {
    # Detect whether run is remote or local
    if (!identical(parent_dir, getwd())) {
        remote <- TRUE
    } else {
        remote <- FALSE
    }

    # Interesting groups, defaults to description
    if (is.null(intgroup)) {
        intgroup <- "description"
    }

    # run_dir
    # Automatic, uses template
    run_dir <- file.path(parent_dir, template)

    # Preliminary run check
    if (!length(dir(run_dir))) {
        stop("run directory failed to load")
    }

    # config_dir
    if (is.null(config_dir)) {
        config_dir <- file.path(run_dir, "config")
    }

    # final_dir
    if (is.null(final_dir)) {
        final_dir <- file.path(run_dir, "final")
    }

    # project_dir
    # Naming scheme is `template/final/YYYY-MM-DD_template`
    project_dir <- file.path(final_dir) %>%
        dir(full.names = TRUE) %>%
        .[grepl(paste0("/\\d{4}-\\d{2}-\\d{2}_[^/]+$"), .)]

    # Sample directories
    sample_dirs <- dir(final_dir, full.names = TRUE)
    names(sample_dirs) <- basename(sample_dirs)
    # Dated `project_dir` is nested here, need to exclude
    sample_dirs <-
        sample_dirs[!grepl("^\\d{4}-\\d{2}-\\d{2}_", names(sample_dirs))]

    # Data versions
    if (file.exists(file.path(project_dir, "data_versions.csv"))) {
        data_versions <- file.path(project_dir, "data_versions.csv") %>%
            read_csv(col_types = cols())
    } else {
        data_versions <- NULL
    }

    # Detect lane split samples
    if (any(grepl("_L\\d+$", dir(final_dir)))) {
        lane_split <- TRUE
    } else {
        lane_split <- FALSE
    }

    run <- list(
        # Input required
        template = template,
        organism = organism,  # can automate by parsing YAML
        # Input recommended
        intgroup = intgroup,
        parent_dir = normalizePath(parent_dir),
        # Automatic for standard runs
        run_dir = normalizePath(run_dir),
        config_dir = normalizePath(config_dir),
        final_dir = normalizePath(final_dir),
        project_dir = normalizePath(project_dir),
        sample_dirs = sample_dirs,
        # Fully automatic
        lane_split = lane_split,
        date = Sys.Date(),
        session_info = sessionInfo(),
        data_versions = data_versions,
        remote = remote,
        wd = getwd()
    )

    check_run(run)

    # Only create new project structure when loading a remote run
    if (isTRUE(remote)) {
        create_new_project()
        save(run, file = "data/run.rda")
    }

    return(run)
}
