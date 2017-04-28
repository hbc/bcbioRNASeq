# Default project naming scheme is `template/final/YYYY-MM-DD_template`

#' Load bcbio run
#'
#' We recommend loading the \code{bcbio-nextgen} run as a remote connection over
#' \code{sshfs}. This requires setting \code{parent_dir} in the function.
#'
#' @author Michael Steinbaugh
#'
#' @param final_dir Path to final output directory. This path is set when
#'   running \code{bcbio_nextgen -w template}.
#' @param intgroup Character vector of interesting groups. First entry is used
#'   for plot colors during quality control (QC) analysis. Entire vector is used
#'   for PCA and heatmap QC functions.
#' @param organism Organism name, following Ensembl/Biomart conventions. Must be
#'   lowercase and one word (e.g. hsapiens). This will be detected automatically
#'   for common reference genomes.
#' @param metadata Optional custom metadata file to import
#'
#' @return \code{bcbio-nextgen} run object
#' @export
load_bcbio_run <- function(
    final_dir = "final",
    intgroup = "description",
    organism = NULL,
    metadata = NULL) {
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

    # Custom metadata
    if (!is.null(metadata)) {
        metadata <- read_metadata(metadata)
    }

    # Detect lane split FASTQ samples
    if (any(str_detect(dir(final_dir), "_L\\d+$"))) {
        lane_split <- TRUE
    } else {
        lane_split <- FALSE
    }

    # YAML summary
    yaml <- read_yaml(file.path(project_dir, "project-summary.yaml"))

    # Detect organism, if not specified
    if (!is.null(yaml)) {
        # Use the genome build of the first sample to match
        organism <- detect_organism(yaml$samples[[1]]$genome_build)
    }
    if (is.null(organism)) {
        stop("Organism is required for metadata queries")
    }

    # Program versions
    programs <- file.path(project_dir, "programs.txt") %>%
        read_delim(",", col_names = c("program", "version"))

    # Data versions
    if (file.exists(file.path(project_dir, "data_versions.csv"))) {
        data_versions <- file.path(project_dir, "data_versions.csv") %>%
            read_csv(col_types = cols())
    } else {
        data_versions <- NULL
    }

    run <- list(
        final_dir = final_dir,
        project_dir = project_dir,
        run_date = as.Date(run_date),
        today_date = Sys.Date(),
        wd = getwd(),
        hpc = detect_hpc(),
        template = template,
        sample_dirs = sample_dirs,
        organism = organism,
        intgroup = intgroup,
        metadata = metadata,
        lane_split = lane_split,
        yaml = yaml,
        programs = programs,
        data_versions = data_versions,
        session_info = sessionInfo()
    )

    check_run(run)
    return(run)
}
