#' Ensembl annotations
#'
#' Download transcript-to-gene annotations from
#' [Ensembl](http://www.ensembl.org/). This function also defines broad classes,
#' which are used in quality control analysis.
#'
#' @author Rory Kirchner
#' @author Michael Steinbaugh
#'
#' @param run bcbio-nextgen run.
#'
#' @return Data frame.
ensembl_annotations <- function(run) {
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
    getBM(mart = ensembl,
          attributes = c("ensembl_gene_id",
                         "ensembl_transcript_id",
                         "external_gene_name",
                         "gene_biotype",
                         "chromosome_name")) %>%
        as_tibble %>%
        group_by(!!sym("ensembl_gene_id")) %>%
        arrange(!!sym("ensembl_transcript_id"), .by_group = TRUE) %>%
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
