# Internal read functions ====
#' Read metadata
#'
#' @keywords internal
#' @author Michael Steinbaugh
#'
#' @param file Metadata file. CSV and XLSX formats are supported.
#' @param pattern Apply grep pattern matching to samples
#' @param pattern_col Column in data frame used for pattern subsetting
#' @param lanes Number of lanes used to split the samples into technical
#'   replicates (`_LXXX`) suffix.
#'
#' @return Metadata data frame.
read_metadata <- function(
    file,
    pattern = NULL,
    pattern_col = "description",
    lanes = NULL) {
    if (!file.exists(file)) {
        stop("File not found")
    }

    if (grepl("\\.xlsx$", file)) {
        metadata <- read_excel(file)
    } else {
        metadata <- read_csv(file)
    }

    # First column must be the FASTQ file name
    names(metadata)[1] <- "file_name"

    metadata <- metadata %>%
        as_tibble %>%
        # Strip all NA rows and columns
        remove_na %>%
        # Make names snake_case
        snake %>%
        # Remove rows with no description
        filter(!is.na(.data$description))

    # Lane split, if desired
    if (is.numeric(lanes)) {
        lane <- paste0("L", str_pad(1:lanes, 3, pad = "0"))
        metadata <- metadata %>%
            group_by(!!sym("file_name")) %>%
            expand_(~lane) %>%
            left_join(metadata, by = "file_name") %>%
            ungroup %>%
            mutate(file_name = paste(.data$file_name, .data$lane, sep = "_"),
                   description = .data$file_name)
    }

    # Subset by pattern, if desired
    if (!is.null(pattern)) {
        metadata <- metadata[str_detect(metadata[[pattern_col]], pattern), ]
    }

    # Convert to data frame, coerce to factors, and set rownames
    metadata %>%
        mutate_all(factor) %>%
        mutate(file_name = as.character(.data$file_name),
               description = as.character(.data$description)) %>%
        arrange(!!sym("description")) %>%
        as.data.frame %>%
        set_rownames(.$description)
}






# Read bcbio-nextgen run output ====
#' Read data from bcbio-nextgen run
#'
#' Read RNA-seq counts using [tximport()]. Currently supports
#' [salmon](https://combine-lab.github.io/salmon/) (**recommended**) and
#' [sailfish](http://www.cs.cmu.edu/~ckingsf/software/sailfish/).
#'
#' @rdname read_bcbio
#' @keywords internal
#'
#' @author Michael Steinbaugh
#' @author Rory Kirchner
#'
#' @param run [bcbioRnaDataSet].
#' @param samples Character vector of sample names in bcbio final directory to
#'   include. If `NULL` (default), all samples will be imported.
#' @param grep Use grep pattern to match samples. This will override the
#'   `samples` parameter if set.
#'
#' @return Counts saved in [tximport] list object.
#' @export
#'
#' @seealso
#' - [tximport::tximport()].
read_bcbio_counts <- function(
    run,
    grep = NULL,
    samples = NULL) {
    # Check for small RNA-seq analysis
    if (run$analysis == "srnaseq") {
        return(list(counts=read_bcbio_file(run, "counts_mirna.tsv", "mirna") %>% as.matrix))
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
    per_sample_dirs <- list.dirs(
        run$sample_dirs[1], full.names = FALSE, recursive = FALSE)
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
#' @param ... Passthrough parameters for [read_csv()] or [read_tsv()].
#'
#' @return Miscellaneous data frame.
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

    # Coerce tibble to data frame
    if (is_tibble(data)) {
        data <- as.data.frame(data)
    }

    # Set row names, if desired
    if (!is.null(row_names)) {
        rownames(data) <- data[[row_names]]
        data[[row_names]] <- NULL
    }

    data %>% remove_na
}



#' @rdname read_bcbio
#' @description Read sample metadata from YAML.
#' @return Metadata data frame.
read_bcbio_metadata <- function(run) {
    read_bcbio_samples_yaml(run, metadata) %>%
        mutate_all(factor) %>%
        mutate(description = as.character(.data$description)) %>%
        arrange(!!sym("description")) %>%
        as.data.frame %>%
        set_rownames(.$description)
}



#' @rdname read_bcbio
#' @description Read summary metrics from YAML. These statistics are only
#'   generated for a standard RNA-seq run with aligned counts. Fast RNA-seq mode
#'   with lightweight counts (pseudocounts) doesn't output the same metrics into
#'   the YAML.
#' @return Summary statistics data frame.
read_bcbio_metrics <- function(run) {
    # [fix] Return without metadata join? May work better downstream?
    yaml <- read_bcbio_samples_yaml(run, summary, metrics)
    if (is.null(yaml)) {
        return(NULL)
    }

    characters <- yaml[, c("description", "quality_format", "sequence_length")]
    numerics <- yaml[, setdiff(colnames(yaml), colnames(characters))] %>%
        mutate_all(as.numeric)
    metrics <- bind_cols(characters, numerics)

    # Return with metadata joined
    if (is.null(run$metadata)) {
        stop("Metadata missing")
    }
    left_join(run$metadata, metrics, by = "description")
}



#' @rdname read_bcbio
#' @description Read bcbio sample information from YAML.
#' @param ... Nested operator keys supplied as dot objects.
#' @return Sample information data frame.
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
        nested <- samples[[a]][[keys]] %>% snake
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

    # List can be coerced to data frame using [data.table::rbindlist()] or
    # [dplyr::bind_rows()]. Some YAML files will cause [bind_rows()] to throw
    # `Column XXX can't be converted from integer to character` errors on
    # numeric data, whereas this doesn't happen with [rbindlist()].
    rbindlist(list) %>%
        as.data.frame %>%
        remove_na %>%
        tidy_select("description", everything()) %>%
        arrange(!!sym("description")) %>%
        set_rownames(.$description)
}






# Small RNA-seq read functions ====
#' Create isomiRs object from bcbio output
#'
#' Read bcbio sample information from YAML to get isomiR object.
#'
#' @rdname read_smallrna
#' @author Lorena Patano
#' @keywords internal
#'
#' @param rna [bcbioRnaDataSet].
read_smallrna_counts <- function(rna) {
    # [fix] Better way to handle sample_dirs than by piping in via metadata?
    run <- metadata(rna)
    fns <- file.path(run$sample_dirs,
                     paste(names(run$sample_dirs),
                           "mirbase-ready.counts",
                           sep = "-"))
    names(fns) <- names(run$sample_dirs)
    message("Reading miRNA count files...")
    bcbio(rna, type = "isomirs") <- IsomirDataSeqFromFiles(
        files = fns[rownames(run$metadata)],
        coldata = run$metadata,
        design = ~run$groups_of_interest)
    rna
}
