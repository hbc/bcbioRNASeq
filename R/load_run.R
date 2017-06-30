#' Load bcbio-nextgen run
#'
#' Simply point to the final upload directory output by
#' [bcbio-nextgen](https://bcbio-nextgen.readthedocs.io/), and this function
#' will take care of the rest. It automatically imports RNA-seq counts,
#' metadata, and program versions used.
#'
#' @note When working in RStudio, we recommend connecting to the bcbio-nextgen
#'   run directory as a remote connection over
#'   [sshfs](https://github.com/osxfuse/osxfuse/wiki/SSHFS).
#'
#' @author Michael Steinbaugh, Lorena Patano
#'
#' @param upload_dir Path to final upload directory. This path is set when
#'   running `bcbio_nextgen -w template`.
#' @param analysis Analysis type. Supports RNA-seq (`rnaseq`; **default**) or
#'   small RNA-seq (`srnaseq`).
#' @param interesting_groups Character vector of interesting groups. First entry
#'   is used for plot colors during quality control (QC) analysis. Entire vector
#'   is used for PCA and heatmap QC functions.
#' @param sample_metadata_file (*Optional*). Custom metadata file containing
#'   sample information. Otherwise defaults to sample metadata saved in the YAML
#'   file.
#' @param experiment_name Experiment name.
#' @param principal_investigator Principal investigator.
#' @param researcher Researcher who performed the experiment.
#' @param email Email for follow-up correspondence.
#'
#' @return [bcbioRNADataSet].
#' @export
load_run <- function(
    upload_dir = "final",
    analysis = "rnaseq",
    interesting_groups = "description",
    design = NULL,
    contrast = NULL,
    sample_metadata_file = NULL) {
    # Directory paths ----
    # Check connection to final upload directory
    if (!dir.exists(upload_dir)) {
        stop("Final upload directory failed to load")
    }
    upload_dir <- normalizePath(upload_dir)
    # Find most recent nested project_dir (normally only 1)
    project_dir_pattern <- "^(\\d{4}-\\d{2}-\\d{2})_([^/]+)$"
    project_dir <- dir(upload_dir,
                       pattern = project_dir_pattern,
                       full.names = FALSE,
                       recursive = FALSE)
    if (length(project_dir) != 1L) {
        stop("Uncertain about project directory location")
    }
    message(project_dir)
    match <- str_match(project_dir, project_dir_pattern)
    run_date <- match[[2L]] %>% as.Date
    template <- match[[3L]]
    project_dir <- file.path(upload_dir, project_dir)


    # Analysis ====
    supported_analyses <- c("rnaseq", "srnaseq")
    if (!analysis %in% supported_analyses) {
        stop(paste("Supported analyses:", toString(supported_analyses)))
    }


    # Project summary YAML ====
    yaml_file <- file.path(project_dir, "project-summary.yaml")
    if (!file.exists(yaml_file)) {
        stop("YAML project summary missing")
    }
    message("Reading project summary YAML")
    yaml <- read_yaml(yaml_file)


    # Sample names ====
    # Obtain the samples (and their directories) from the YAML
    sample_names <- vapply(
        yaml[["samples"]],
        function(x) x[["description"]],
        character(1L)) %>% sort
    if (!identical(basename(sample_dirs), sample_names)) {
        stop("Sample name assignment mismatch")
    }
    sample_dirs <- file.path(upload_dir, sample_names) %>%
        set_names(sample_names)
    message(paste(length(sample_dirs), "samples detected"))


    # Genome ====
    # Use the genome build of the first sample to match
    genome_build <- yaml[["samples"]][[1L]][["genome_build"]]
    organism <- .detect_organism(genome_build)
    message(paste("Genome:", organism, genome_build))
    annotable <- .annotable(genome_build)
    tx2gene <- .tx2gene(project_dir, genome_build)


    # Sequencing lanes ====
    lane_pattern <- "_L(\\d{3})"
    if (any(str_detect(sample_dirs, lane_pattern))) {
        lanes <- str_match(names(sample_dirs), lane_pattern) %>%
            .[, 2L] %>%
            unique %>%
            length
        message(paste(
            lanes, "sequencing lane detected", "(technical replicates)"))
    } else {
        # TODO Check that downstream functions don't use NULL
        lanes <- 1L
    }


    # User-defined custom metadata ====
    sample_metadata <- .sample_metadata_file(sample_metadata_file)


    # Sample metrics ====
    # Note that sample metrics used for QC plots are not currently generated
    # when using fast RNA-seq workflow. This depends upon MultiQC and aligned
    # counts generated with STAR.
    metrics <- .yaml_metrics(yaml)


    # Data versions and programs ====
    data_versions <- .data_versions(project_dir)
    programs <- .programs(project_dir)


    # Log files ====
    bcbio_nextgen <- read_lines(
        file.path(project_dir, "bcbio-nextgen.log"))
    bcbio_nextgen_commands <- read_lines(
        file.path(project_dir, "bcbio-nextgen-commands.log"))


    # Metadata ====
    metadata <- SimpleList(
        analysis = analysis,
        upload_dir = upload_dir,
        sample_dirs = sample_dirs,
        project_dir = project_dir,
        template = template,
        run_date = run_date,
        interesting_groups = interesting_groups,
        organism = organism,
        genome_build = genome_build,
        ensembl_version = ensembl_version,
        annotable = annotable,
        tx2gene = tx2gene,
        lanes = lanes,
        yaml_file = yaml_file,
        yaml = yaml,
        metrics = metrics,
        sample_metadata_file = sample_metadata_file,
        sample_metadata = sample_metadata,
        data_versions = data_versions,
        programs = programs,
        bcbio_nextgen = bcbio_nextgen,
        bcbio_nextgen_commands = bcbio_nextgen_commands,
        experiment_name = experiment_name,
        principal_investigator = principal_investigator,
        researcher = researcher,
        email = email)


    # tximport ====
    tximport <- .tximport(sample_dirs, tx2gene = tx2gene)


    # SummarizedExperiment ====
    if (!is.null(sample_metadata)) {
        col_data <- sample_metadata
    } else {
        col_data <- .yaml_metadata(yaml)
    }
    se <- .summarized_experiment(
        tximport = tximport,
        col_data = col_data,
        row_data = annotable,
        metadata = metadata)


    # bcbioRNADataSet ====
    bcb <- new("bcbioRNADataSet", se)
    bcbio(bcb, "tximport") <- tximport
    bcb
}
