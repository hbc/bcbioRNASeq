#' Read metadata
#'
#' @param file Metadata file. CSV and XLSX formats are supported.
#' @param pattern Apply grep pattern matching to samples
#' @param pattern_col Column in data frame used for pattern subsetting
#' @param lanes Number of lanes used to split the samples into technical
#'   replicates. This assumes the file names don't include `_L001` suffix.
#'   Therefore, it renames both the file name and description columns to match
#'   the bcbio YAML. This workflow is commonly used by Harvard Biopolymers
#'   Facility. We can decide to either combine counts at the server level using
#'   `cat` in bash, or we can run DESeq2 later by pooling the counts with
#'   `deseq_lane_pool()`. We may want to deprecate this method in the future and
#'   simply combine counts at the server level for all lane split runs.
#' @return Metadata data frame
#' @export
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
        set_names_snake %>%
        .[!is.na(.$file_name), ] %>%
        .[!is.na(.$description), ] %>%
        .[order(.$description), ]

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

    # Convert to data frame and set rownames
    metadata <- metadata %>%
        as.data.frame %>%
        set_rownames(.$description)

    return(metadata)
}
