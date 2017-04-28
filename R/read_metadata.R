#' Read metadata
#'
#' @param csv Comma separated values file
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
    csv,
    pattern = NULL,
    pattern_col = "description",
    lanes = NULL) {
    metadata <- read_csv(csv, col_types = cols()) %>%
        set_names_snake %>%
        arrange(!!sym("description"))
    # First column must be the FASTQ file name
    names(metadata)[1] <- "file_name"

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

    return(metadata)
}
