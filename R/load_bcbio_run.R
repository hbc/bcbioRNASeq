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
#' @param parent_dir Parent directory containing the run folder
#' @param config_dir Set the config output directory, if non-standard
#' @param final_dir Set the final output directory, if non-standard
#' @param factor Condition factor vector
#' @param intgroup Internal grouping character (string) used for plot colors.
#'   Defaults to the first entry in the factor vector if not specified.
#'
#' @return bcbio list object with directory paths
#' @export
#'
#' @examples
#' \dontrun{
#' load_bcbio_run("illumina_rnaseq")
#' }
load_bcbio_run <-
    function(template,
             organism,
             parent_dir = NULL,
             config_dir = NULL,
             final_dir = NULL,
             factor,
             intgroup = NULL) {
        if (is.null(parent_dir)) {
            parent_dir <- getwd()
        }

        # `bcbio-nextgen` run
        run_dir <- file.path(parent_dir, template)
        if (!length(dir(run_dir))) {
            stop("Run directory failed to load.")
        }

        if (is.null(config_dir)) {
            config_dir <- file.path(run_dir, "config")
        }
        if (!length(dir(config_dir))) {
            stop("Config directory failed to load.")
        }

        if (is.null(final_dir)) {
            final_dir <- file.path(run_dir, "final")
        }
        if (!length(dir(final_dir))) {
            stop("Final directory failed to load.")
        }

        # Project naming scheme is `template/final/YYYY-MM-DD_template`
        project_dir <- file.path(final_dir) %>%
            dir(full.names = TRUE) %>%
            .[grepl(paste0("/\\d{4}-\\d{2}-\\d{2}_[^/]+$"), .)]
        if (!length(dir(project_dir))) {
            stop("Project directory failed to load.")
        }

        # Work on documenting this in better detail
        if (is.null(factor)) {
            stop("A condition factor vector is required.")
        }

        if (is.null(intgroup)) {
            intgroup <- factor[1]
        }

        # Detect lane split samples
        if (any(grepl("_L\\d+$", dir(final_dir)))) {
            lane_split <- TRUE
        } else {
            lane_split <- FALSE
        }

        # Create directories for R project
        dir.create("data", showWarnings = FALSE)
        dir.create("results", showWarnings = FALSE)

        bcbio <- list(
            parent_dir = normalizePath(parent_dir),
            run_dir = normalizePath(run_dir),
            config_dir = normalizePath(config_dir),
            final_dir = normalizePath(final_dir),
            project_dir = normalizePath(project_dir),
            organism = organism,
            factor = factor,
            intgroup = intgroup,
            lane_split = lane_split
        )

        save(bcbio, file = "data/bcbio.rda")
        return(bcbio)
    }
