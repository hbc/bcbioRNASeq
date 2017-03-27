#' Run DESeq2 with pooled technical replicates
#'
#' @author Michael Steinbaugh
#'
#' @import DESeq2
#'
#' @param txi tximport list object
#' @param intgroup intgroup
#' @param design Design formula
#' @return DESeqDataSet object
run_deseq_pooled <- function(
    txi,
    intgroup,
    design,
    metadata) {
    counts <- as.matrix(txi$counts)
    stem <- gsub("_L\\d+$", "", colnames(counts)) %>% unique %>% sort

    counts <- lapply(seq_along(stem), function(a) {
        grep <- paste0("^", stem[a], "_L\\d+$")
        counts %>% .[, grepl(grep, colnames(.))] %>% rowSums
    }) %>%
        set_names(stem) %>%
        do.call(cbind, .)

    metadata <- import_metadata(bcbio,
                                intgroup = intgroup,
                                lane_split = FALSE)

    dds <- DESeq2::DESeqDataSetFromMatrix(
        counts,
        colData = metadata,
        design = design) %>%
        DeSeq2::DESeq(.)

    return(dds)
}
