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
.ensembl <- function(genome_build) {
    # Remap genome build versions to annotables
    if (genome_build == "hg19") {
        genome_build <- "grch37"
    } else if (genome_build == "hg38") {
        genome_build <- "grch38"
    } else if (genome_build == "mm10") {
        genome_build <- "grcm38"
    }

    annotables <- .annotables(genome_build)

    list <- SimpleList()
    list[["gene"]] <- annotables[["gene"]] %>%
        tidy_select(c("ensgene", "symbol", "chr", "biotype", "description")) %>%
        distinct %>%
        arrange(!!sym("ensgene")) %>%
        as.data.frame %>%
        column_to_rownames("ensgene") %>%
        DataFrame
    list[["tx2gene"]] <- annotables[["tx2gene"]] %>%
        arrange(!!sym("enstxp")) %>%
        as.data.frame

    list
}



#' @rdname ensembl
.annotables <- function(genome_build) {
    envir <- as.environment("package:annotables")
    gene <- get(genome_build, envir = envir)
    tx2gene <- paste(genome_build, "tx2gene", sep = "_") %>%
        get(envir = envir)
    list(gene = gene, tx2gene = tx2gene)
}
