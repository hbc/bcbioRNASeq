#' Read data from bcbio-nextgen run
#'
#' @rdname read_bcbio
#'
#' @author Michael Steinbaugh
#' @author Rory Kirchner
#'
#' @param run \code{bcbio-nextgen} run



#' @rdname read_bcbio
#' @description Read RNA-Seq counts using \code{tximport}
#'
#' @param type The type of software used to generate the abundances. Follows the
#'   conventions of \code{tximport}.
#' @param grep Apply grep pattern matching to samples
#' @param samples Specify the names of samples in bcbio final directory to
#'   input. If \code{NULL} (default), all samples will be imported.
#'
#' @return txi \code{tximport} list
#' @export
read_bcbio_counts <- function(
    run,
    type = "salmon",
    grep = NULL,
    samples = NULL) {
    check_run(run)
    if (!type %in% c("salmon", "sailfish")) {
        stop("Unsupported counts input format")
    }

    # Sample name grep pattern matching. Run above `samples` to override.
    if (!is.null(grep)) {
        samples <- str_subset(names(run$sample_dirs), pattern = grep)
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
        stop(paste(type, "Count files do not exist"))
    }

    # Begin loading of selected counts
    message(paste("loading", type, "counts"))
    message(paste(names(sample_files), collapse = "\n"))

    tx2gene <- read_bcbio_file(run, file = "tx2gene.csv", col_names = FALSE)

    # Import the counts
    # https://goo.gl/h6fm15
    # countsFromAbundance = c("no", "scaledTPM", "lengthScaledTPM")
    # Use `lengthScaledTPM` for `salmon` and `sailfish`
    txi <- tximport(files = sample_files,
                    type = type,
                    tx2gene = tx2gene,
                    importer = read_tsv,
                    countsFromAbundance = "lengthScaledTPM")

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
    read_bcbio_samples_yaml(run, keys = "metadata")
}



#' @rdname read_bcbio
#' @description Read summary metrics from YAML
#' @return Summary statistics data frame
#' @export
read_bcbio_metrics <- function(run) {
    meta <- read_bcbio_metadata(run)
    metrics <- read_bcbio_samples_yaml(run, keys = c("summary", "metrics"))
    left_join(meta, metrics, by = "description")
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
    check_run(run)
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

        # Set empty and NULL values as NA
        empty_as_na <- function(a) {
            if (is.null(a)) {
                a <- NA
            }
            a <- gsub("^$", NA, a)
            return(a)
        }
        nested <- lapply(nested, empty_as_na)
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
