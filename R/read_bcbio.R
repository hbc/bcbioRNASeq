#' Read data from bcbio-nextgen run.
#'
#' @rdname read_bcbio
#'
#' @author Michael Steinbaugh
#' @author Rory Kirchner
#'
#' @param run bcbio-nextgen run.



#' @rdname read_bcbio
#' @description Read RNA-seq counts using \code{\link[tximport]{tximport}}.
#'   Currently supports \href{https://combine-lab.github.io/salmon/}{salmon} and
#'   \href{http://www.cs.cmu.edu/~ckingsf/software/sailfish/}{sailfish}.
#'
#' @param samples Character vector of sample names in bcbio final directory to
#'   include. If \code{NULL} (default), all samples will be imported.
#' @param grep Use grep pattern to match samples. This will override the
#'   \code{samples} parameter if set.
#'
#' @return txi \code{\link[tximport]{tximport}} list.
#' @export
read_bcbio_counts <- function(
    run,
    grep = NULL,
    samples = NULL) {
    # Check for small RNA-seq analysis
    if (run$analysis == "srnaseq") {
        read_bcbio_file(run, "counts_mirna.tsv", "mirna") %>% as.matrix
    }

    # Select the samples to import
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
    # Ensure the sample names are sorted alphabetically
    samples <- sort(samples)

    # Check for count output format, by using the first sample directory
    per_sample_dirs <- list.dirs(run$sample_dirs[1],
                                 full.names = FALSE,
                                 recursive = FALSE)
    if ("salmon" %in% per_sample_dirs) {
        type <- "salmon"
    } else if ("sailfish" %in% per_sample_dirs) {
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
    names(sample_files) %>% toString %>% message

    # Import the counts
    # https://goo.gl/h6fm15
    # countsFromAbundance = c("no", "scaledTPM", "lengthScaledTPM")
    if (type %in% c("salmon", "sailfish")) {
        # Use `lengthScaledTPM` for salmon and sailfish
        countsFromAbundance <- "lengthScaledTPM"
    } else {
        countsFromAbundance <- "no"
    }

    tximport(files = sample_files,
             type = type,
             tx2gene = run$tx2gene,
             importer = read_tsv,
             countsFromAbundance = countsFromAbundance)
}



#' @rdname read_bcbio
#' @description Read a project data file.
#'
#' @param file File name.
#' @param row_names Column identifier to use for row names.
#' @param ... Passthrough parameters for the \code{readr} importers called
#'   internally (\code{\link[readr]{read_csv}}, \code{\link[readr]{read_tsv}}).
#'
#' @return bcbio run data.
#' @export
#'
#' @seealso
#' \code{\link[readr]{read_csv}}, \code{\link[readr]{read_tsv}})
read_bcbio_file <- function(
    run,
    file,
    row_names = NULL,
    ...) {
    # Check that file exists
    filepath <- file.path(run$project_dir, file)
    if (!file.exists(filepath)) {
        return(NULL)
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
        data <- read_csv(filepath, ...)
    } else if (ext == "tsv") {
        data <- read_tsv(filepath, ...)
    } else if (ext == "counts") {
        data <- read_tsv(filepath, ...) %>% as.matrix
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

    data
}



#' @rdname read_bcbio
#' @description Read sample metadata from YAML.
#' @return Metadata data frame.
#' @export
read_bcbio_metadata <- function(run) {
    read_bcbio_samples_yaml(run, metadata)
}



#' @rdname read_bcbio
#' @description Read summary metrics from YAML.
#' @return Summary statistics data frame.
#' @export
read_bcbio_metrics <- function(run) {
    metadata <- run$metadata
    # These statistics are only generated for a standard RNA-seq run with
    # aligned counts. Fast RNA-seq mode with lightweight counts (pseudocounts)
    # doesn't output the same metrics into the YAML.
    metrics <- read_bcbio_samples_yaml(run, summary, metrics)
    left_join(metadata, metrics, by = "description")
}



#' @rdname read_bcbio
#' @description Read bcbio sample information from YAML.
#' @keywords internal
#'
#' @param ... Nested operator keys supplied as dot objects.
#'
#' @export
read_bcbio_samples_yaml <- function(run, ...) {
    keys <- as.character(substitute(list(...)))[-1L]
    yaml <- run$yaml
    if (is.null(yaml)) {
        stop("Run YAML summary is required")
    }
    samples <- yaml$samples
    if (!length(samples)) {
        stop("No sample information in YAML")
    }

    # Check presence of nested keys, otherwise return NULL
    # [fix] Improve the recursion method using sapply in a future update
    if (!keys[1] %in% names(samples[[1]])) {
        return(NULL)
    }
    if (length(keys) > 1) {
        if (!keys[2] %in% names(samples[[1]][[keys[1]]])) {
            return(NULL)
        }
    }

    list <- lapply(seq_along(samples), function(a) {
        nested <- samples[[a]][[keys]] %>%
            # Sanitize names in snake_case
            set_names_snake
        # Set the description
        nested$description <- samples[[a]]$description
        # Remove legacy duplicate `name` identifier
        nested$name <- NULL
        # Correct batch and phenotype YAML
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
        nested
    })

    # [data.table::rbindlist()] coerces integers better than
    # [dplyr::bind_rows()]. Some YAML files will cause [bind_rows()] to throw
    # "Column `appy_severity` can't be converted from integer to character"
    # errors on numeric data.
    rbindlist(list) %>%
        as.data.frame %>%
        # Put description first and sort other colnames alphabetically
        .[order(.$description), ] %>%
        .[, c("description", sort(setdiff(names(.), "description")))] %>%
        set_rownames(.$description)
}



#' Create isomiRs object from bcbio output
#'
#' Read bcbio output to create isomiRs object
#'
#' @rdname read_smallrna_counts
#' @description Read bcbio sample information from YAML to get
#' isomiR object
#'
#' @param bcbiods bcbioRnaDataSet object
#'
#' @export
read_smallrna_counts <- function(bcbiods){
    run <- metadata(bcbiods)
    fns <- file.path(run$sample_dirs,
                     paste(names(run$sample_dirs),
                           "mirbase-ready.counts",
                           sep = "-"))
    names(fns) <- names(run$sample_dirs)
    message("Reading miRNA count files...")
    IsomirDataSeqFromFiles(files = fns[rownames(run$metadata)],
                           coldata = run$metadata,
                           design = ~run$intgroup)
}
