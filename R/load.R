#' Load bcbio-nextgen run
#'
#' Simply point to the final upload directory output by
#' [bcbio-nextgen](https://bcbio-nextgen.readthedocs.io/), and this function
#' will take care of the rest. It automatically imports RNA-seq counts,
#' metadata, and program versions used.
#'
#' When working in RStudio, we recommend connecting to the bcbio-nextgen run
#' directory as a remote connection over
#' [sshfs](https://github.com/osxfuse/osxfuse/wiki/SSHFS).
#'
#' @author Michael Steinbaugh
#' @author Lorena Patano
#'
#' @param upload_dir Path to final upload directory. This path is set when
#'   running `bcbio_nextgen -w template`.
#' @param analysis Analysis type (e.g. `rnaseq` or `srnaseq`).
#' @param intgroup Character vector of interesting groups. First entry is used
#'   for plot colors during quality control (QC) analysis. Entire vector is used
#'   for PCA and heatmap QC functions.
#' @param metadata_file Custom metadata file to import. Otherwise defaults to
#'   sample metadata saved in the YAML file.
#' @param organism Organism name, following Ensembl/Biomart conventions. Must be
#'   lowercase and one word (e.g. hsapiens). This will be detected automatically
#'   for common reference genomes and normally does not need to be set.
#'
#' @return bcbio-nextgen run object (as list), containing counts and metadata.
#' @export
load_run <- function(
    upload_dir = "final",
    analysis = "rnaseq",
    intgroup = "description",
    metadata_file = NULL,
    organism = NULL) {
    # Check connection to upload_dir
    if (!dir.exists(upload_dir)) {
        stop("Final upload directory failed to load")
    }
    upload_dir <- normalizePath(upload_dir)

    # Check analysis type
    if (!analysis %in% c("rnaseq", "srnaseq")) {
        stop("Unsupported analysis type")
    }

    # Find most recent nested project_dir (normally only 1)
    pattern <- "^(\\d{4}-\\d{2}-\\d{2})_([^/]+)$"
    project_dir <- dir(upload_dir,
                       pattern = pattern,
                       full.names = FALSE,
                       recursive = FALSE) %>%
        sort %>% rev %>% .[[1]]
    if (!length(project_dir)) {
        stop("Project directory not found")
    }
    message(project_dir)
    match <- str_match(project_dir, pattern)
    # run_date <- match[2]
    template <- match[3]
    project_dir <- file.path(upload_dir, project_dir)

    # Program versions
    message("Reading program versions...")
    programs <- file.path(project_dir, "programs.txt") %>%
        read_delim(",", col_names = c("program", "version"))

    # Data versions
    if (file.exists(file.path(project_dir, "data_versions.csv"))) {
        message("Reading data versions...")
        data_versions <- read_csv(file.path(project_dir, "data_versions.csv"))
    } else {
        data_versions <- NULL
    }

    # Find the summary YAML file automatically
    yaml_file <- dir(project_dir, pattern = "*.yaml$", full.names = TRUE)
    if (length(yaml_file) != 1) {
        stop("Unsure which YAML file to use")
    }

    # Load the YAML summary file
    message(paste("YAML:", basename(yaml_file)))
    yaml <- read_yaml(yaml_file)

    # Organism class (used with Ensembl)
    if (is.null(organism)) {
        # Use the genome build of the first sample to match
        genome_build <- yaml$samples[[1]]$genome_build
        organism <- detect_organism(genome_build)
    }
    message(paste("Genome:", genome_build))
    message(paste("Organism:", organism))

    # Obtain the samples (and their directories) from the YAML
    sample_names <- sapply(
        yaml$samples, function(x) { x$description }) %>% sort
    sample_dirs <- file.path(upload_dir, sample_names)
    names(sample_dirs) <- sample_names
    message(paste(length(sample_dirs), "samples detected in run"))

    # Detect number of sample lanes, determine if split
    lane_pattern <- "_L(\\d{3})"
    if (any(str_detect(sample_dirs, lane_pattern))) {
        # [fix] message("Sample lanes detected")
        lanes <- str_match(names(sample_dirs), lane_pattern) %>%
            .[, 2] %>% unique %>% length
        message(paste(lanes, "lane replicates per sample detected"))
    } else {
        lanes <- NULL
    }

    run <- list(
        yaml_file = yaml_file,
        yaml = yaml,
        analysis = analysis,
        upload_dir = upload_dir,
        project_dir = project_dir,
        run_date = as.Date(yaml$date),
        today_date = Sys.Date(),
        template = template,
        wd = getwd(),
        hpc = detect_hpc(),
        sample_dirs = sample_dirs,
        intgroup = intgroup,
        metadata_file = metadata_file,
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
        run$metadata <- read_metadata(metadata_file, lanes = lanes)
    } else {
        run$metadata <- read_bcbio_metadata(run)
    }

    # Subset the sample_dirs by the metadata data frame
    run$sample_dirs <- run$sample_dirs[run$metadata$description]
    message(paste(length(run$sample_dirs), "samples matched by metadata"))

    if (analysis == "rnaseq") {
        # Save Ensembl annotations
        message("Downloading Ensembl annotations...")
        run$ensembl <- ensembl(run)
        run$ensembl_version <- listMarts() %>% .[1, 2]

        # Use the tx2gene file output by `bcbio-nextgen`. Alternatively,
        # we can handle this directly in R using the Ensembl annotations
        # obtained with biomaRt instead, if the file doesn't exist. This
        # is currently the case with the fast RNA-seq pipeline, for example.
        run$tx2gene <- read_bcbio_file(
            run, file = "tx2gene.csv", col_names = FALSE)
        if (is.null(run$tx2gene)) {
            run$tx2gene <- tx2gene(run)
        }
    }

    if (analysis %in% c("rnaseq", "srnaseq")) {
        # Read counts using [tximport()]
        run$txi <- read_bcbio_counts(run)
    } else {
        # Placeholder for future workflows
    }

    check_run(run)
    create_project_dirs()
    download_rmarkdown_files()
    run
}
