#' Run DESeq2 with pooled technical replicates
#'
#' @author Michael Steinbaugh
#'
#' @import DESeq2
#'
#' @param txi tximport list object
#' @param intgroup intgroup
#' @param design Design formula
#'
#' @return DESeqDataSet object
#' @export
run_deseq_pooled <- function(
    txi,
    intgroup,
    design) {
    counts <- as.matrix(txi$counts)
    stem <- gsub("_L\\d+$", "", colnames(counts)) %>% unique %>% sort
    
    pooled_counts <- lapply(seq_along(stem), function(a) {
        grep <- paste0("^", stem[a], "_L\\d+$")
        counts %>% .[, grepl(grep, colnames(.))] %>% rowSums
    }) %>%
        set_names(stem) %>%
        do.call(cbind, .) %>%
        # Counts must be integers, otherwise DESeq will error
        round
    
    metadata <- import_metadata(bcbio,
                                intgroup = intgroup,
                                lane_split = FALSE)
    
    dds <- DESeq2::DESeqDataSetFromMatrix(
        pooled_counts,
        colData = metadata,
        design = design
    ) %>% DESeq2::DESeq(.)
    
    return(dds)
}
