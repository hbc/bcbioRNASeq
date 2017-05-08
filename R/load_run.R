#' Load bcbio run
#'
#' We recommend loading the \code{bcbio-nextgen} run as a remote connection over
#' \code{sshfs}. This requires setting \code{parent_dir} in the function.
#'
#' @author Michael Steinbaugh
#'
#' @param yaml_file YAML file saved during final upload from
#'   \code{bcbio-nextgen}
#' @param analysis Analysis type (e.g. \code{rnaseq} or \code{srnaseq})
#' @param intgroup Character vector of interesting groups. First entry is used
#'   for plot colors during quality control (QC) analysis. Entire vector is used
#'   for PCA and heatmap QC functions.
#' @param metadata_file Custom metadata file to import. Otherwise defaults to
#'   sample metadata saved in the YAML file.
#' @param organism Organism name, following Ensembl/Biomart conventions. Must be
#'   lowercase and one word (e.g. hsapiens). This will be detected automatically
#'   for common reference genomes and normally does not need to be set.
#'
#' @return \code{bcbio-nextgen} run object
#' @export
load_run <- function(
    yaml_file = "project-summary.yaml",
    analysis = "rnaseq",
    intgroup = "description",
    metadata_file = NULL,
    organism = NULL) {
    # Load the YAML summary file
    if (file.exists(yaml_file)) {
        yaml_file <- normalizePath(yaml_file)
        message(paste("YAML:", yaml_file))
        yaml <- read_yaml(yaml_file)
    } else {
        stop("Run YAML file not found")
    }

    # Recursve two levels up from YAML to get final upload directory base
    project_dir <- yaml_file %>% dirname
    upload_dir <- yaml_file %>% dirname %>% dirname
    message(paste(dir(upload_dir), collapse = "\n"))

    # Obtain the samples (and their directories) from the YAML
    sample_names <- sapply(yaml$samples,
                            function(x) { x$description }) %>% sort
    sample_dirs <- file.path(upload_dir, sample_names)
    names(sample_dirs) <- sample_names

    # Organism class (used with Ensembl)
    if (is.null(organism)) {
        # Use the genome build of the first sample to match
        organism <- detect_organism(yaml$samples[[1]]$genome_build)
    }

    # Detect number of sample lanes, determine if split
    lane_pattern <- "_L(\\d{3})"
    if (any(str_detect(sample_dirs, lane_pattern))) {
        message("Sample lanes detected")
        lanes <- str_match(names(sample_dirs), lane_pattern) %>%
            .[, 2] %>% unique %>% length
    } else {
        lanes <- NULL
    }

    # Program versions
    message("Reading program versions...")
    programs <- file.path(project_dir, "programs.txt") %>%
        read_delim(",", col_names = c("program", "version"))

    # Data versions
    if (file.exists(file.path(project_dir, "data_versions.csv"))) {
        message("Reading data versions...")
        data_versions <- file.path(project_dir, "data_versions.csv") %>%
            read_csv(col_types = cols())
    } else {
        data_versions <- NULL
    }

    run <- list(
        yaml_file = yaml_file,
        yaml = yaml,
        analysis = analysis,
        upload_dir = upload_dir,
        project_dir = project_dir,
        run_date = as.Date(yaml$date),
        today_date = Sys.Date(),
        wd = getwd(),
        hpc = detect_hpc(),
        sample_dirs = sample_dirs,
        intgroup = intgroup,
        metadata = NULL,
        organism = organism,
        ensembl = NULL,
        ensembl_version = NULL,
        lanes = lanes,
        programs = programs,
        data_versions = data_versions,
        session_info = sessionInfo())

    # Store metadata
    if (!is.null(metadata_file)) {
        # Read from custom data frame
        run$metadata <- read_metadata(metadata_file, lanes = lanes)
    } else{
        # Read from YAML
        run$metadata <- read_bcbio_samples_yaml(run, keys = "metadata")
    }

    # Subset the sample_dirs by the metadata data frame
    run$sample_dirs <- run$sample_dirs[run$metadata$description]

    # Save Ensembl annotations
    if (analysis == "rnaseq") {
        run$ensembl <- ensembl_annotations(
            run,
            attributes = c("external_gene_name",
                           "description",
                           "gene_biotype"))
        run$ensembl_version <- listMarts() %>% .[1, 2]
    }

    check_run(run)
    return(run)
}
