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
#' @param groups_of_interest Character vector of interesting groups. First entry
#'   is used for plot colors during quality control (QC) analysis. Entire vector
#'   is used for PCA and heatmap QC functions.
#' @param custom_metadata_file (*Optional*). Custom metadata file containing
#'   sample information. Otherwise defaults to sample metadata saved in the YAML
#'   file.
#'
#' @return [bcbioRnaDataSet].
#' @export
load_run <- function(
    upload_dir = "final",
    analysis = "rnaseq",
    groups_of_interest = "description",
    custom_metadata_file = NULL) {
    # Directory paths ----
    # Check connection to final upload directory
    if (!dir.exists(upload_dir)) {
        stop("Final upload directory failed to load")
    }
    upload_dir <- normalizePath(upload_dir)
    # Find most recent nested project_dir (normally only 1)
    upload_dir_pattern <- "^(\\d{4}-\\d{2}-\\d{2})_([^/]+)$"
    project_dir <- dir(upload_dir,
                       pattern = upload_dir_pattern,
                       full.names = FALSE,
                       recursive = FALSE) %>%
        sort %>% rev %>% .[[1]]
    if (!length(project_dir)) {
        stop("Project directory missing")
    }
    message(project_dir)
    match <- str_match(project_dir, upload_dir_pattern)
    run_date <- match[2] %>% as.Date
    template <- match[3]
    project_dir <- file.path(upload_dir, project_dir)


    # Analysis ----
    supported_analyses <- c("rnaseq", "srnaseq")
    if (!analysis %in% supported_analyses) {
        stop(paste("Supported analyses:", toString(supported_analyses)))
    }


    # Project summary YAML ----
    yaml_file <- file.path(project_dir, "project-summary.yaml")
    if (!file.exists(yaml_file)) {
        stop("YAML project summary missing")
    }
    message("Reading project summary YAML")
    yaml <- read_yaml(yaml_file)


    # Sample names ----
    # Obtain the samples (and their directories) from the YAML
    sample_names <- vapply(
        yaml[["samples"]],
        function(x) x[["description"]],
        character(1L)) %>% sort
    sample_dirs <- file.path(upload_dir, sample_names) %>%
        set_names(sample_names)
    if (!identical(basename(sample_dirs), sample_names)) {
        stop("Sample name assignment mismatch")
    }
    message(paste(length(sample_dirs), "samples detected"))


    # Genome ----
    # Use the genome build of the first sample to match
    genome_build <- yaml[["samples"]][[1L]][["genome_build"]]
    organism <- .detect_organism(genome_build)
    message(paste("Genome:", organism, genome_build))
    # Ensembl annotations from annotables
    ensembl <- .ensembl(genome_build)


    # Sequencing lanes ----
    lane_pattern <- "_L(\\d{3})"
    if (any(str_detect(sample_dirs, lane_pattern))) {
        lanes <- str_match(names(sample_dirs), lane_pattern) %>%
            .[, 2L] %>% unique %>% length
        message(paste(lanes, "lane replicates per sample detected"))
    } else {
        # [TODO] Check that downstream functions don't use NULL
        lanes <- 1L
    }


    # SummarizedExperiment ----
    # [TODO] Pass in custom metadata
    custom_metadata <- DataFrame()
    se <- .SummarizedExperiment(
        txi = .tximport(sample_dirs, tx2gene = ensembl[["tx2gene"]]),
        colData = .yaml_metadata(yaml),
        rowData = ensembl[["gene"]],
        metadata = custom_metadata)


    # bcbioRnaDataSet ----
    bcb <- new(
        "bcbioRnaDataSet", se,
        analysis = analysis,
        groups_of_interest = groups_of_interest,
        # design
        # contrast
        upload_dir = upload_dir,
        project_dir = project_dir,
        sample_dirs = sample_dirs,
        run_date = run_date,
        load_date = Sys.Date(),
        template = template,
        organism = organism,
        genome_build = genome_build,
        yaml_file = yaml_file,
        yaml = yaml,
        custom_metadata_file = custom_metadata_file,
        lanes = lanes,  # [TODO] Set lanes to 1 instead of null
        data_versions = .data_versions(project_dir),
        program_versions = .program_versions(project_dir),
        wd = getwd(),
        hpc = detect_hpc(),
        session_info = sessionInfo())


    # Optional slots ----
    if (!is.null(alt_counts)) {
        # [TODO] Update to slot alt_counts
        # bcb@alt_counts <- alt_counts
    } else if (!is.null(custom_metadata_file)) {
        bcb@custom_metadata_file <- custom_metadata_file
        # [TODO] Add workflow to process custom metadata
        custom_metadata <- DataFrame()
        bcb@custom_metadata <- custom_metadata
    }
}
