#' Detect the organism from the genome build name
#'
#' @keywords internal
#' @author Michael Steinbaugh
#'
#' @param genome_build Genome build
#'
#' @return Organism string
#' @export
#'
#' @examples
#' detect_organism("hg19")
#' detect_organism("mm10")
#' detect_organism("WBcel235")
#' detect_organism("BDGP6")
detect_organism <- function(genome_build) {
    if (str_detect(genome_build, "^(hg|GRCh)\\d+")) {
        "hsapiens"
    } else if (str_detect(genome_build, "^mm\\d+")) {
        "mmusculus"
    } else if (str_detect(genome_build, "^WBcel\\d+")) {
        "celegans"
    } else if (str_detect(genome_build, "^BDGP\\d+")) {
        "dmelanogaster"
    } else if (str_detect(genome_build, "^Zv\\d+")) {
        "drerio"
    } else if (str_detect(genome_build, "^ASM\\d+")) {
        "spombe"
    } else if (str_detect(genome_build, "^(MB|MG)\\d+")) {
        "ecoli"
    } else {
        stop("Failed to detect organism")
    }
}



#' Read metadata
#'
#' @keywords internal
#' @author Michael Steinbaugh
#'
#' @param file Metadata file. CSV and XLSX formats are supported.
#' @param pattern Apply grep pattern matching to samples
#' @param pattern_col Column in data frame used for pattern subsetting
#' @param lanes Number of lanes used to split the samples into technical
#'   replicates (\code{_LXXX}) suffix.
#'
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
    metadata %>%
        as.data.frame %>%
        set_rownames(.$description)
}



res_contrast_name <- function(res) {
    mcols(res)[2, 2] %>%
        str_replace("^.*:\\s", "")
}



tx2gene <- function(run) {
    run$ensembl %>%
        as.data.frame %>%
        .[, c("ensembl_transcript_id", "ensembl_gene_id")] %>%
        set_rownames(.$ensembl_transcript_id)
}



tx2gene_annotations <- function(run) {
    run$ensembl %>%
        ungroup %>%
        mutate(ensembl_transcript_id = NULL) %>%
        distinct
}
