# aggregate_replicates ====
#' @rdname aggregate_replicates
#'
#' @param object Object.
#' @param pattern Grep pattern to match lane identifiers in file name.
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
    stem <- str_replace(colnames(object), pattern, "") %>% unique %>% sort
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
#' @export
#' @seealso A new [DESeqDataSet] is returned using [DESeqDataSetFromMatrix()].
setMethod("aggregate_replicates", "DESeqDataSet", function(
    object,
    pattern = "_L\\d+") {
    # Get the stored design formula
    design <- design(dds)

    # Pool the lane split technical replicates
    # `countData` is defined in [DESeqDataSetFromMatrix()]
    message("Aggregating raw counts slotted in DESeqDataSet")
    countData <- counts(dds, normalized = FALSE) %>%
        aggregate_replicates

    # Mutate the colData metadata to pooled samples
    colData <- colData(dds) %>%
        as.data.frame %>%
        mutate(file_name = str_replace(.data$file_name, pattern, ""),
               description = str_replace(.data$description, pattern, ""),
               lane = NULL,
               # [DESeq()]-generated columns
               # replaceable = NULL,
               sizeFactor = NULL) %>%
        distinct %>%
        set_rownames(.$description)

    # Check that the new colData matches the counts matrix
    if (!identical(colData$file_name, colnames(countData))) {
        stop("File name mismatch in colData and countData")
    }

    # Replace file names in pooled count matrix with description
    colnames(countData) <- colData$description
    if (!identical(rownames(colData), colnames(countData))) {
        stop("Description mismatch in colData and countData")
    }

    message("Reloading DESeqDataSet using DESeqDataSetFromMatrix")
    DESeqDataSetFromMatrix(
        countData = countData,
        colData = colData,
        design = design) %>% DESeq
})



# bcbio ====
#' @rdname bcbio
#' @export
setMethod("bcbio", "bcbioRnaDataSet", function(object, type = "counts") {
    if (!type %in% c("counts", "abundance", "length", "alt_counts")) {
        stop("Unsupported type")
    }
    assays(object)[[type]]
})

#' @rdname bcbio
#' @exportMethod "bcbio<-"
setReplaceMethod(
    "bcbio",
    signature(object = "bcbioRnaDataSet", value = "matrix"),  # ANY
    function(object, type = "counts", value) {
        if (!type %in% c("counts", "abundance", "length", "alt_counts")) {
            stop("Unsupported type")
        }
        assays(object)[[type]] <- value
        validObject(object)
        object
    })



# melt_log10 ====
.melt_log10 <- function(counts) {
    counts %>%
        as.data.frame %>%
        rownames_to_column %>%
        melt(id = 1) %>%
        set_names(c("rowname",  # ensembl_gene_id
                    "colname",  # description
                    "counts")) %>%
        # Filter zero counts
        filter(.data[["counts"]] > 0) %>%
        # log10 transform
        mutate(counts = log10(.data[["counts"]]),
               # [melt()] sets colnames as factor
               colname = as.character(.data[["colname"]]))
}

# [TODO] Add TMM count matrix support
#' @rdname melt_log10
#' @export
setMethod("melt_log10", "bcbioRnaDataSet", function(
    object,
    normalized = FALSE,
    interesting_groups = NULL) {
    if (isTRUE(normalized)) {
        counts <- tpm(object)
    } else {
        counts <- raw_counts(object)
    }

    if (is.null(interesting_groups)) {
        interesting_groups <- metadata(object)[["interesting_groups"]]
    }

    metadata <- colData(object) %>%
        as.data.frame %>%
        # Use "colname" here, corresponding to colData above
        rownames_to_column("colname") %>%
        tidy_select(!!!syms(c("colname", interesting_groups)))

    if (!identical(colnames(counts), metadata[["colname"]])) {
        stop("Sample description mismatch between counts and metadata")
    }

    .melt_log10(counts) %>%
        left_join(metadata, by = "colname") %>%
        rename(ensgene = .data[["rowname"]],
               description = .data[["colname"]])
})


#' @rdname melt_log10
#' @export
setMethod("melt_log10", "DESeqDataSet", function(
    object,
    normalized = FALSE,
    interesting_groups = NULL) {
    if (isTRUE(normalized)) {
        counts <- DESeq2::counts(object, normalized = TRUE)
    } else {
        counts <- DESeq2::counts(object, normalized = FALSE)
    }

    metadata <- colData(dds)
    if (!is.null(interesting_groups)) {
    }

})
# DESeqDataSet(normalized = TRUE)
# DESeqTransform



# metadata_table ====
#' @rdname metadata_table
#' @export
setMethod("metadata_table", "bcbioRnaDataSet", function(object, ...) {
    object %>%
        colData %>%
        as.data.frame %>%
        remove_rownames %>%
        kable(caption = "Sample metadata", ...)
})

# DESeqDataSet



# metrics ====
#' @rdname metrics
#' @export
setMethod("metrics", "bcbioRnaDataSet", function(object) {
    metrics <- metadata(bcb)[["metrics"]]
    if (is.null(metrics)) return(NULL)
    cbind(colData(bcb), metrics) %>% as.data.frame
})



# raw_counts ====
#' @rdname raw_counts
#' @export
setMethod("raw_counts", "bcbioRnaDataSet", function(object) {
    assays(object)[["counts"]]
})

# DESeqDataSet



# sample_dirs ====
#' @rdname sample_dirs
#' @export
setMethod("sample_dirs", "bcbioRnaDataSet", function(object) {
    metadata(object)[["sample_dirs"]]
})



# tpm ====
#' @rdname tpm
#' @export
setMethod("tpm", "bcbioRnaDataSet", function(object) {
    assays(object)[["abundance"]]
})
