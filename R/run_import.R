#' Import data from bcbio-nextgen run
#'
#' @rdname run_import
#'
#' @author Michael Steinbaugh
#' @author Rory Kirchner
#'
#' @param run \code{bcbio-nextgen}
#' @param save Save data



#' @rdname run_import
#' @description Import metadata
#'
#' @param grep Apply grep pattern matching to samples
#' @param pool Pool lane split samples
#'
#' @return Metadata data frame
#' @export
import_metadata <- function(
    run,
    grep = NULL,
    save = FALSE,
    pool = FALSE) {
    check_run(run)
    metadata <- list.files(run$config_dir,
                           pattern = ".csv$",
                           full.names = TRUE) %>%
        read_csv(col_types = cols()) %>%
        set_names_snake

    if (isTRUE(run$lane_split) & !isTRUE(pool)) {
        # Lane splitting. This assumes the YAML descriptions won't match the
        # `_L00[1-4]` suffix. Therefore, it renames both the samplename and
        # description columns to match the bcbio server output. This workflow is
        # used by Harvard Biopolymers Facility. We can decide to either combine
        # counts at the server level using `cat` in bash, or we can run DESeq2
        # later by pooling the counts with `deseq_lane_pool()`. We may want to
        # deprecate this method in the future and simply combine counts at the
        # server level for all lane split runs.
        lane <- paste0("L", stringr::str_pad(1:4, 3, pad = "0"))

        metadata <- metadata %>%
            group_by_(.dots = "samplename") %>%
            expand_(~lane) %>%
            left_join(metadata, by = "samplename") %>%
            ungroup
        metadata$samplename <- paste(metadata$samplename, lane, sep = "_")

        # Check against the sample_dirs
        description_match <- identical(metadata$description,
                                       names(run$sample_dirs))
        samplename_match <- identical(metadata$samplename,
                                      names(run$sample_dirs))

        # Replace description for lanesplit samples
        if (!isTRUE(description_match) & isTRUE(samplename_match)) {
            metadata$description <- metadata$samplename
        }
    } else if (isTRUE(pool)) {
        metadata$description <- metadata$samplename
    }

    if (!is.null(grep)) {
        metadata <- metadata[grepl(grep, metadata$description), ]
    }

    metadata <- metadata %>%
        arrange_(.dots = "description") %>%
        as.data.frame %>%
        set_rownames(.$description)

    if (isTRUE(save)) {
        save(metadata, file = "data/metadata.rda")
        write_csv(metadata, "meta/metadata.csv")
    }

    return(metadata)
}



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



#' @rdname run_import
#' @description Import RNA-Seq counts using \code{tximport}
#'
#' @param type The type of software used to generate the abundances. Follows the
#'   conventions of \code{tximport}.
#' @param samples Specify the names of samples in bcbio final directory to
#'   input. If \code{NULL} (default), all samples will be imported.
#'
#' @return txi \code{tximport} list
#' @export
import_counts <- function(
    run,
    type = "salmon",
    grep = NULL,
    samples = NULL) {
    check_run(run)
    if (!type %in% c("salmon", "sailfish")) {
        stop("unsupported counts input format")
    }

    # Sample name grep pattern matching. Run above `samples` to override.
    if (!is.null(grep)) {
        samples <- stringr::str_subset(
            names(run$sample_dirs), pattern = grep)
        if (!length(samples)) {
            stop("grep string didn't match any samples")
        }
    }

    if (is.null(samples)) {
        samples <- names(run$sample_dirs)
    } else {
        samples <- sort(unique(na.omit(samples)))
    }

    # Support for salmon and sailfish file structure
    sample_files <- file.path(run$sample_dirs[samples],
                              type, "quant", "quant.sf")
    names(sample_files) <- names(run$sample_dirs[samples])
    if (!all(file.exists(sample_files))) {
        stop(paste(type, "count files do not exist"))
    }

    # Begin loading of selected counts
    message(paste("loading", type, "counts"))
    message(paste(names(sample_files), collapse = "\n"))
    # writeLines(sample_files)

    tx2gene <- import_file(run, file = "tx2gene.csv", col_names = FALSE)

    # Import the counts
    # https://goo.gl/h6fm15
    # countsFromAbundance = c("no", "scaledTPM", "lengthScaledTPM")
    # Use `lengthScaledTPM` for `salmon` and `sailfish`
    txi <- tximport(files = sample_files,
                    type = type,
                    tx2gene = tx2gene,
                    reader = read_tsv,
                    countsFromAbundance = "lengthScaledTPM")

    return(txi)
}



#' @rdname run_import
#' @description Import a project data file
#'
#' @param file File name
#' @param row_names Column identifier to use for row names
#' @param ... Optional parameters for \code{readr}
#'
#' @return bcbio run data
#' @export
import_file <- function(
    run,
    file,
    row_names = NULL,
    ...) {
    check_run(run)

    # Check that file exists
    filepath <- file.path(run$project_dir, file)
    if (!file.exists(filepath)) {
        stop("file could not be found")
    }

    # Detect file extension
    if (grepl("\\.[a-z]+$", file)) {
        # ext <- gsub("^.*\\.([a-z]+)$", "\\1", file)
        ext <- str_match(file, "\\.([a-z]+)$")[2]
    } else {
        stop("file does not have an extension")
    }

    # File import
    if (ext == "csv") {
        data <- read_csv(filepath, col_types = cols(), ...)
    } else if (ext == "tsv") {
        data <- read_tsv(filepath, col_types = cols(), ...)
    } else if (ext == "counts") {
        data <- read_tsv(filepath, col_types = cols(), ...) %>% as.matrix
    } else {
        stop("unsupported file extension")
    }

    # Coerce tibble to data frame. Might want to disable this in a future update
    # as tibble becomes more standardized and the rownames issue is sorted out.
    if (is_tibble(data)) {
        data <- as.data.frame(data)
    }

    # Set row names
    if (!is.null(row_names)) {
        rownames(data) <- data[[row_names]]
        data[[row_names]] <- NULL
    }

    return(data)
}
