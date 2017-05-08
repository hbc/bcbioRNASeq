#' Read data from bcbio-nextgen run
#'
#' @rdname read_bcbio
#'
#' @author Michael Steinbaugh
#' @author Rory Kirchner
#'
#' @param run \code{bcbio-nextgen} run



#' @rdname read_bcbio
#' @description Read RNA-seq counts using \code{tximport}. Currently supports
#'   \href{https://combine-lab.github.io/salmon/}{salmon} and
#'   \href{http://www.cs.cmu.edu/~ckingsf/software/sailfish/}{sailfish}.
#'
#' @param samples Character vector of sample names in bcbio final directory to
#'   include. If \code{NULL} (default), all samples will be imported.
#' @param grep Use grep pattern to match samples. This will override the
#'   \code{samples} parameter if set.
#'
#' @return txi \code{tximport} list
#' @export
read_bcbio_counts <- function(
    run,
    grep = NULL,
    samples = NULL) {
    check_run(run)
    if (!is.null(grep)) {
        # grep pattern matching against sample names
        samples <- str_subset(names(run$sample_dirs), pattern = grep)
        if (!length(samples)) {
            stop("grep pattern didn't match any samples")
        }
    } else if (!is.null(samples)) {
        samples <- sort(unique(na.omit(samples)))
    } else {
        # Include all samples by default
        samples <- names(run$sample_dirs)
    }

    # Check for count output format, by using the first sample directory
    per_sample_dirs <- list.dirs(run$sample_dirs[1], recursive = FALSE)
    if (any(str_detect(per_sample_dirs, "salmon"))) {
        type <- "salmon"
    } else if (any(str_detect(per_sample_dirs, "sailfish"))) {
        type <- "sailfish"
    } else {
        stop("Unsupported counts output format")
    }

    # Locate quant.sf file for salmon or sailfish output
    if (type %in% c("salmon", "sailfish")) {
        sample_files <- list.files(
            file.path(run$sample_dirs[samples], type),
            pattern = "quant.sf",
            full.names = TRUE,
            recursive = TRUE)
        # Check that count files exist for all samples
        if (length(sample_files) != length(samples)) {
            stop(paste(type, "counts not found for all samples"))
        }
    }

    # Assign descriptions to sample files
    names(sample_files) <- names(run$sample_dirs[samples])

    # Begin loading of selected counts
    message(paste("Reading", type, "counts..."))
    message(paste(names(sample_files), collapse = "\n"))

    # Use the tx2gene file output by `bcbio-nextgen`. Alternatively,
    # we could handle this directly in R using biomaRt instead.
    tx2gene <- read_bcbio_file(run, file = "tx2gene.csv", col_names = FALSE)

    # Import the counts
    # https://goo.gl/h6fm15
    # countsFromAbundance = c("no", "scaledTPM", "lengthScaledTPM")
    if (type %in% c("salmon", "sailfish")) {
        # Use `lengthScaledTPM` for `salmon` and `sailfish`
        countsFromAbundance <- "lengthScaledTPM"
    } else {
        # default
        countsFromAbundance <- "no"
    }
    txi <- tximport(files = sample_files,
                    type = type,
                    tx2gene = tx2gene,
                    importer = read_tsv,
                    countsFromAbundance = countsFromAbundance)

    return(txi)
}



#' @rdname read_bcbio
#' @description Read a project data file
#'
#' @param file File name
#' @param row_names Column identifier to use for row names
#' @param ... Optional parameters for \code{readr}
#'
#' @return bcbio run data
#' @export
read_bcbio_file <- function(
    run,
    file,
    row_names = NULL,
    ...) {
    check_run(run)

    # Check that file exists
    filepath <- file.path(run$project_dir, file)
    if (!file.exists(filepath)) {
        stop("File could not be found")
    }

    # Detect file extension
    if (grepl("\\.[a-z]+$", file)) {
        # ext <- gsub("^.*\\.([a-z]+)$", "\\1", file)
        ext <- str_match(file, "\\.([a-z]+)$")[2]
    } else {
        stop("File does not have an extension")
    }

    # File import
    if (ext == "csv") {
        data <- read_csv(filepath, col_types = cols(), ...)
    } else if (ext == "tsv") {
        data <- read_tsv(filepath, col_types = cols(), ...)
    } else if (ext == "counts") {
        data <- read_tsv(filepath, col_types = cols(), ...) %>% as.matrix
    } else {
        stop("Unsupported file extension")
    }

    # Set row names, if desired
    if (!is.null(row_names)) {
        # Coerce tibble to data frame. Setting rownames on a tibble is now
        # deprecated. If rownames are essential, object must be a data frame.
        if (is_tibble(data)) {
            data <- as.data.frame(data)
        }
        rownames(data) <- data[[row_names]]
        data[[row_names]] <- NULL
    }

    return(data)
}



#' @rdname read_bcbio
#' @description Read sample metadata from YAML
#' @return Metadata data frame
#' @export
read_bcbio_metadata <- function(run) {
    check_run(run)
    read_bcbio_samples_yaml(run, keys = "metadata")
}



#' @rdname read_bcbio
#' @description Read summary metrics from YAML
#' @return Summary statistics data frame
#' @export
read_bcbio_metrics <- function(run) {
    check_run(run)
    metadata <- run$metadata
    metrics <- read_bcbio_samples_yaml(run, keys = c("summary", "metrics"))
    left_join(metadata, metrics, by = "description")
}



#' @rdname read_bcbio
#' @description Read bcbio sample information from YAML
#' @keywords internal
#'
#' @param keys Nested operator keys that should be supplied as an ordered
#'     character vector, recursing a level down for each entry
#'
#' @export
read_bcbio_samples_yaml <- function(run, keys) {
    # Don't run integrity checks here, used to save into the run object
    yaml <- run$yaml
    if (is.null(yaml)) {
        stop("Run YAML summary is required")
    }
    samples <- yaml$samples
    if (!length(samples)) {
        stop("No sample information in YAML")
    }

    list <- lapply(seq_along(samples), function(a) {
        nested <- samples[[a]][[keys]] %>%
            # Sanitize names in snake_case
            set_names_snake
        # Set the description
        nested$description <- samples[[a]]$description
        # Remove legacy duplicate `name` identifier
        nested$name <- NULL

        # Correct batch and phenotype for metadata, if selected
        if (rev(keys)[1] == "metadata") {
            # Fix empty batch and phenotype
            if (is.null(nested$batch)) {
                nested$batch <- NA
            }
            if (length(nested$phenotype)) {
                if (grepl("^$", nested$phenotype)) {
                    nested$phenotype <- NA
                }
            }
        }

        # Coerce numerics for metrics, if selected
        if (rev(keys)[1] == "metrics") {
            char_vec <- c("description",
                          "quality_format",
                          "sequence_length")
            characters <- nested[names(nested) %in% char_vec]
            numerics <- setdiff(names(nested),
                                names(characters)) %>%
                sort %>%
                nested[.] %>%
                lapply(as.numeric)
            nested <- append(characters, numerics)
        }
        return(nested)
    })

    df <- bind_rows(list) %>%
        as.data.frame %>%
        .[order(.$description), ] %>%
        # Put description first and sort everything else
        .[, c("description", sort(setdiff(names(.), "description")))] %>%
        set_rownames(.$description)
    return(df)
}
