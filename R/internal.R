#' Detect the organism from the genome build name
#'
#' @keywords internal
#' @author Michael Steinbaugh
#'
#' @param genome_build Genome build
#'
#' @return Organism string
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



#' Return metadata intgroups as factor
#'
#' @author Michael Steinbaugh
#' @keywords internal
#'
#' @param run bcbio-nextgen run.
#'
#' @export
intgroup_as_factor <- function(run) {
    import_tidy_verbs()
    run$metadata %>%
        select(!!!syms(c("description", run$intgroup))) %>%
        mutate_all(factor) %>%
        as.data.frame %>%
        set_rownames(.$description) %>%
        .[, run$intgroup]
}



#' Query Ensembl
#'
#' @keywords internal
#' @author Michael Steinbaugh
#' @author Rory Kirchner
#'
#' @param run bcbio run.
#'
#' @return Tibble grouped by `ensembl_gene_id` with nested
#'   `ensembl_transcript_id`.
query_ensembl <- function(run) {
    import_tidy_verbs()
    # Broad class definitions
    coding <- c("protein_coding")
    decaying <- c("non_stop_decay", "nonsense_mediated_decay")
    noncoding <- c("known_ncrna", "lincRNA", "non_coding")
    srna <- c("miRNA", "misc_RNA", "ribozyme", "rRNA", "scaRNA", "scRNA",
              "snoRNA", "snRNA", "sRNA")
    ensembl <- useEnsembl(
        biomart = "ensembl",
        dataset = paste(run$organism, "gene_ensembl", sep = "_"))
    nested_tx2gene <- getBM(
        mart = ensembl,
        attributes = c("ensembl_gene_id",
                       "ensembl_transcript_id")) %>%
        as_tibble %>%
        group_by(!!sym("ensembl_gene_id")) %>%
        arrange(!!sym("ensembl_transcript_id"), .by_group = TRUE) %>%
        nest_(key_col = "ensembl_transcript_id",
              nest_cols = "ensembl_transcript_id")
    gene_level <- getBM(
        mart = ensembl,
        attributes = c("ensembl_gene_id",
                       "external_gene_name",
                       "description",
                       "gene_biotype",
                       "chromosome_name"))
    left_join(nested_tx2gene, gene_level, by = "ensembl_gene_id") %>%
        mutate(broad_class = case_when(
            tolower(.data$chromosome_name) == "mt" ~ "mito",
            # Fix to match Drosophila genome (non-standard)
            grepl("mito", .data$chromosome_name) ~ "mito",
            grepl("pseudo", .data$gene_biotype) ~ "pseudo",
            grepl("TR_", .data$gene_biotype) ~ "TCR",
            grepl("IG_", .data$gene_biotype) ~ "IG",
            .data$gene_biotype %in% srna ~ "small",
            .data$gene_biotype %in% decaying ~ "decaying",
            .data$gene_biotype %in% noncoding ~ "noncoding",
            .data$gene_biotype %in% coding ~ "coding",
            TRUE ~ "other"))
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
        snake %>%
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
