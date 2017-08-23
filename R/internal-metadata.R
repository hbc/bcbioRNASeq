.interestingColData <- function(object) {
    colData(object) %>%
        as.data.frame %>%
        .[, c("sampleName", .interestingGroups(object)), drop = FALSE]
}



.interestingGroup <- function(object) {
    metadata(object) %>%
        .[["interestingGroups"]] %>%
        .[[1L]] %>%
        as.name
}



.interestingGroups <- function(object) {
    metadata(object)[["interestingGroups"]]
}



.metaPriorityCols <- function(meta) {
    meta %>%
        as("tibble") %>%
        # Sanitize `sampleID` into valid names
        mutate(sampleID = make.names(.data[["sampleName"]])) %>%
        .[, unique(c(metaPriorityCols, sort(colnames(.))))] %>%
        arrange(!!!syms(metaPriorityCols))
}



.metaFactors <- function(meta) {
    meta %>%
        mutate_if(!colnames(.) %in% metaPriorityCols, factor)
}



.uniqueMetrics <- function(object) {
    dropCols <- c(metaPriorityCols, "name")
    metrics <- metadata(object)[["metrics"]] %>%
        as.data.frame %>%
        set_rownames(.[["sampleID"]]) %>%
        .[, setdiff(colnames(.), dropCols), drop = FALSE]
    # Find metrics columns with unique values
    keepCols <- lapply(colnames(metrics), function(a) {
        if (length(unique(metrics[, a])) > 1L) a
    }) %>%
        unlist %>%
        .[!is.null(.)]
    metrics[, keepCols, drop = FALSE]
}
