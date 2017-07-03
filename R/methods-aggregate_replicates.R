#' Aggregate lane-split technical replicates
#'
#' Frequently RNA-seq experiments are performed with technical replicates
#' split across flow cell lanes. This generic facilitates quick aggregation
#' of counts across the flow cells.
#'
#' @rdname aggregate_replicates
#' @docType methods
#'
#' @author Michael Steinbaugh
#'
#' @param object Object.
#' @param pattern Grep pattern to match lane identifiers in file name.
#' @param ... Additional arguments.
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
        mutate(file_name = str_replace(.data[["file_name"]], pattern, ""),
               description = str_replace(.data[["description"]], pattern, ""),
               lane = NULL,
               sizeFactor = NULL) %>%
        distinct %>%
        set_rownames(.[["description"]])

    # Check that the new col_data matches the counts matrix
    if (!identical(col_data[["file_name"]], colnames(count_data))) {
        stop("File name mismatch in col_data and count_data")
    }

    # Replace file names in pooled count matrix with description
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
