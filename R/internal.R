# General ====
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



#' Get contrast name from [DESeqResults]
#'
#' @keywords internal
#' @author Michael Steinbaugh
#'
#' @param res [DESeqResults].
#'
#' @return Contrast name string.
res_contrast_name <- function(res) {
    mcols(res)[2, 2] %>%
        str_replace("^.*:\\s", "")
}



#' [DESeq] [colData]
#'
#' Selects `intgroup`-defined interesting groups automatically.
#'
#' @keywords internal
#' @author Michael Steinbaugh
#'
#' @param dds [DESeqDataSet] or [DESeqTransform].
#' @param intgroup Interesting groups.
#'
#' @return [colData()] data frame.
select_intgroup_coldata <- function(dds, intgroup) {
    colData(dds) %>%
        as.data.frame %>%
        tidy_select(!!quo(intgroup))
}






# Integrity checks ====
#' Perform object integrity checks
#'
#' @rdname integrity_checks
#' @keywords internal
#' @author Michael Steinbaugh
#'
#' @param run [bcbioRnaDataSet].
#' @param stop Stop upon check failure.
#'
#' @export
check_run <- function(run) {
    if (!is.list(run)) {
        stop("Run object is not a list")
    }
    if (is.null(run$yaml)) {
        stop("Run does contain YAML info, please resave")
    }
    if (is.null(run$metadata)) {
        stop("Run does not contain metadata, please resave")
    }
    if (run$analysis == "rnaseq" & is.null(run$ensembl)) {
        message("Run does not contain Ensembl annotations")
    }

    # Metadata format and description name checks
    if (is_tibble(run$metadata) | !is.data.frame(run$metadata)) {
        stop("Run metadata not saved as a data frame")
    }
    if (!identical(rownames(run$metadata), run$metadata$description)) {
        stop("Run metadata rownames must match the description")
    }
    if (!identical(names(run$sample_dirs), run$metadata$description)) {
        stop("Run metadata descriptions don't match sample directories")
    }
    if (!any(sapply(run$metadata, is.factor))) {
        stop("Run metadata does not contain factors")
    }
}



#' @rdname integrity_checks
#' @param dds [DESeqDataSet].
check_dds <- function(dds, stop = TRUE) {
    if (class(dds)[1] == "DESeqDataSet") {
        TRUE
    } else {
        if (isTRUE(stop)) {
            stop("DESeqDataSet required")
        } else {
            FALSE
        }
    }
}



#' @rdname integrity_checks
#' @param dt [DESeqTransform].
check_dt <- function(dt, stop = TRUE) {
    if (class(dt)[1] == "DESeqTransform") {
        TRUE
    } else {
        if (isTRUE(stop)) {
            stop("DESeqDataSet required")
        } else {
            FALSE
        }
    }
}



#' @rdname integrity_checks
#' @param res [DESeqResults].
check_res <- function(res, stop = TRUE) {
    if (class(res)[1] == "DESeqResults") {
        TRUE
    } else {
        if (isTRUE(stop)) {
            stop("DESeqResults required")
        } else {
            FALSE
        }
    }
}







# Ensembl ====
#' Ensembl annotations
#'
#' Transcript-level annotations are downloaded from
#' [Ensembl](http://www.ensembl.org/) and saved as a tibbled grouped by
#' `ensembl_gene_id` with nested `ensembl_transcript_id` in the bcbio run
#' object. We also define broad classes, which are used in quality control
#' analysis.
#'
#' @rdname ensembl
#' @keywords internal
#'
#' @author Michael Steinbaugh
#' @author Rory Kirchner
#'
#' @param run [bcbioRnaDataSet].
#'
#' @return Tibble grouped by `ensembl_gene_id` with nested
#'   `ensembl_transcript_id`.
ensembl <- function(run) {
    # Broad class definitions
    if ( run$genome_build %in% c("GRCh37", "hg19", "hg38", "mm10"))
        return(.annotable(run$genome_build))
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
    .merge_annotation(nested_tx2gene, gene_level)
}


.annotable <- function(build){
    # [! Fix] use string version for distinct
    if (build %in% c("GRCh37", "hg19")){
        gene_table <- annotables::grch37 %>% distinct(ensgene,.keep_all = TRUE)
        tx_table <- annotables::grch37_tx2gene
    }else if (build == "hg38"){
        gene_table <- annotables::grch38 %>% distinct(ensgene,.keep_all = TRUE)
        tx_table <- annotables::grch38_tx2gene
    }else if (build == "mm10"){
        gene_table <- annotables::grcm38 %>% distinct(ensgene,.keep_all = TRUE)
        tx_table <- annotables::grcm38_tx2gene
    }else{
        stop("annotation not supported")
    }
    names(tx_table) <- c("ensembl_transcript_id", "ensembl_gene_id")
    tx_table <- tx_table %>% as_tibble %>%
        group_by(!!sym("ensembl_gene_id")) %>%
        arrange(!!sym("ensembl_transcript_id"), .by_group = TRUE) %>%
        nest_(key_col = "ensembl_transcript_id",
              nest_cols = "ensembl_transcript_id")
    names(gene_table) <- c("ensembl_gene_id",
                           "entrez_id",
                           "external_gene_name",
                           "chromosome_name",
                           "start", "end", "strand",
                           "gene_biotype",
                           "description")
    .merge_annotation(tx_table, gene_table)
}

# support ensembl
.merge_annotation <- function(nested_tx2gene, gene_level){
    coding <- c("protein_coding")
    decaying <- c("non_stop_decay", "nonsense_mediated_decay")
    noncoding <- c("known_ncrna", "lincRNA", "non_coding")
    srna <- c("miRNA", "misc_RNA", "ribozyme", "rRNA", "scaRNA", "scRNA",
              "snoRNA", "snRNA", "sRNA")

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


#' @rdname ensembl
#' @description Gene-level annotations with `description` and `gene_biotype`,
#'   and `broad_class`.
#' @return Tibble arranged by `ensembl_gene_id`.
#' @export
gene_level_annotations <- function(run) {
    run$ensembl %>%
        mutate(ensembl_transcript_id = NULL) %>%
        ungroup %>%
        distinct %>%
        arrange(!!sym("ensembl_gene_id")) %>%
        as.data.frame %>%
        set_rownames(.$ensembl_gene_id)
}



#' @rdname ensembl
#' @description Transcript-to-gene annotations, to be used with [tximport()].
#' @return Tibble with `ensembl_transcript_id` and `ensembl_gene_id`.
#' @export
tx2gene <- function(run) {
    run$ensembl %>%
        unnest_("ensembl_transcript_id") %>%
        tidy_select(!!!syms(c("ensembl_transcript_id", "ensembl_gene_id"))) %>%
        ungroup %>%
        arrange(!!sym("ensembl_transcript_id")) %>%
        as.data.frame %>%
        set_rownames(.$ensembl_transcript_id)
}
