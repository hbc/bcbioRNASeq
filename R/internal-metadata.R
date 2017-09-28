.interestingColData <- function(object) {
    cols <- unique(c("sampleName", .interestingGroups(object)))
    colData(object) %>%
        as.data.frame() %>%
        .[, cols, drop = FALSE]
}



.interestingGroup <- function(object) {
    metadata(object) %>%
        .[["interestingGroups"]] %>%
        .[[1L]] %>%
        as.character()
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
