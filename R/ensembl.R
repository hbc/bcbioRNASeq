#' Ensembl annotations
#'
#' @details Transcript-level annotations are downloaded from
#'   [Ensembl](http://www.ensembl.org/) and saved as a tibbled grouped by
#'   `ensembl_gene_id` with nested `ensembl_transcript_id` in the bcbio run
#'   object. We also define broad classes, which are used in quality control
#'   analysis.
#'
#' @rdname ensembl
#'
#' @author Michael Steinbaugh
#'
#' @param run bcbio-nextgen run.



#' @rdname ensembl
#' @description Gene-level annotations with `description` and `gene_biotype`,
#'   and `broad_class`.
#' @return Tibble arranged by `ensembl_gene_id`.
#' @export
gene_level_annotations <- function(run) {
    import_tidy_verbs()
    run$ensembl %>%
        mutate(ensembl_transcript_id = NULL) %>%
        ungroup %>%
        distinct %>%
        arrange(!!sym("ensembl_gene_id"))
}



#' @rdname ensembl
#' @description Transcript-to-gene annotations, to be used with [tximport()].
#' @return Tibble with `ensembl_transcript_id` and `ensembl_gene_id`.
#' @export
tx2gene <- function(run) {
    import_tidy_verbs()
    run$ensembl %>%
        unnest_("ensembl_transcript_id") %>%
        select(!!!syms(c("ensembl_transcript_id", "ensembl_gene_id"))) %>%
        ungroup %>%
        arrange(!!sym("ensembl_transcript_id"))
}
