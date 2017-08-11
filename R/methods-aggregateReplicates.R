#' Aggregate Lane-Split Technical Replicates
#'
#' Frequently RNA-seq experiments are performed with technical replicates
#' split across flow cell lanes. This generic facilitates quick aggregation
#' of counts across the flow cells.
#'
#' @rdname aggregateReplicates
#' @name aggregateReplicates
#'
#' @param pattern Grep pattern to match lane identifiers in sample name.
#'
#' @return Object of same class, with pooled technical replicates.
NULL



# Methods ====
#' @rdname aggregateReplicates
#' @export
setMethod("aggregateReplicates", "matrix", function(
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
        setNames(stem) %>%
        do.call(cbind, .) %>%
        # [round()] here otherwise [DESeq()] will fail on this matrix
        round
})



#' @rdname aggregateReplicates
#' @note [DESeqDataSet] is returned using [DESeqDataSetFromMatrix()].
#' @export
setMethod("aggregateReplicates", "DESeqDataSet", function(
    object,
    pattern = "_L\\d+") {
    # Get the stored design formula
    design <- design(object)

    # Pool the lane split technical replicates
    message("Aggregating raw counts slotted in DESeqDataSet")
    countData <- counts(object, normalized = FALSE) %>%
        aggregateReplicates

    # Mutate the colData metadata to pooled samples
    colData <- colData(object) %>%
        as.data.frame %>%
        mutate(sampleID = str_replace(.data[["sampleID"]], pattern, ""),
               sampleName = str_replace(.data[["sampleName"]], pattern, ""),
               lane = NULL,
               sizeFactor = NULL) %>%
        distinct %>%
        set_rownames(.[["sampleID"]])

    # Check that the new colData matches the counts matrix
    if (!identical(colData[["sampleName"]], colnames(countData))) {
        stop("Sample name mismatch in colData and countData")
    }

    # Replace sample names in pooled count matrix with description
    colnames(countData) <- colData[["description"]]
    if (!identical(rownames(colData), colnames(countData))) {
        stop("Description mismatch in colData and countData")
    }

    message("Reloading DESeqDataSet using DESeqDataSetFromMatrix")
    DESeqDataSetFromMatrix(
        countData = countData,
        colData = colData,
        design = design) %>% DESeq
})
