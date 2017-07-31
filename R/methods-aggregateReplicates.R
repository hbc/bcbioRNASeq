#' Aggregate Lane-Split Technical Replicates
#'
#' Frequently RNA-seq experiments are performed with technical replicates
#' split across flow cell lanes. This generic facilitates quick aggregation
#' of counts across the flow cells.
#'
#' @rdname aggregate_replicates
#' @author Michael Steinbaugh
#'
#' @param pattern Grep pattern to match lane identifiers in sample name.
#'
#' @return Object of same class, with pooled technical replicates.
#' @export
setMethod("aggregate_replicates", "matrix", function(
    object,
    pattern = "_L\\d+") {
    # Obtain the unique pooled sample names
    if (!all(grepl(pattern, colnames(object)))) {
        stop("Lane pattern didn't match all samples")
    }
    stem <- str_replace(colnames(object), pattern, "") %>%
        unique %>%
        sort
    # Perform [rowSums()] on the matching columns per sample
    lapply(seq_along(stem), function(a) {
        object %>%
            .[, grepl(paste0("^", stem[a], pattern), colnames(.))] %>%
            rowSums
    }) %>%
        set_names(stem) %>%
        do.call(cbind, .) %>%
        # [round()] here otherwise [DESeq()] will fail on this matrix
        round
})



#' @rdname aggregate_replicates
#' @note [DESeqDataSet] is returned using [DESeqDataSetFromMatrix()].
#' @export
setMethod("aggregate_replicates", "DESeqDataSet", function(
    object,
    pattern = "_L\\d+") {
    # Get the stored design formula
    design <- design(object)

    # Pool the lane split technical replicates
    message("Aggregating raw counts slotted in DESeqDataSet")
    count_data <- counts(object, normalized = FALSE) %>%
        aggregate_replicates

    # Mutate the colData metadata to pooled samples
    col_data <- colData(object) %>%
        as.data.frame %>%
        mutate(sample_id = str_replace(.data[["sample_id"]], pattern, ""),
               sample_name = str_replace(.data[["sample_name"]], pattern, ""),
               lane = NULL,
               sizeFactor = NULL) %>%
        distinct %>%
        set_rownames(.[["sample_id"]])

    # Check that the new col_data matches the counts matrix
    if (!identical(col_data[["sample_name"]], colnames(count_data))) {
        stop("Sample name mismatch in col_data and count_data")
    }

    # Replace sample names in pooled count matrix with description
    colnames(count_data) <- col_data[["description"]]
    if (!identical(rownames(col_data), colnames(count_data))) {
        stop("Description mismatch in col_data and count_data")
    }

    message("Reloading DESeqDataSet using DESeqDataSetFromMatrix")
    DESeqDataSetFromMatrix(
        countData = count_data,
        colData = col_data,
        design = design) %>% DESeq
})
