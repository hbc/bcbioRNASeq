# TODO Add support for loading aligned counts instead of salmon

#' @importFrom magrittr set_colnames
#' @importFrom readr read_tsv
.featureCounts <- function(projectDir) {
    inform("Reading STAR featureCounts aligned counts")
    file <- path(projectDir, "combined.counts")
    assert_all_are_existing_files(file)
    read_tsv(featureCountsFile) %>%
        as.data.frame() %>%
        # Sanitize colnames into valid names
        set_colnames(
            gsub(
                x = make.names(colnames(.), unique = TRUE),
                pattern = "\\.",
                replacement = "_")
        ) %>%
        column_to_rownames("id") %>%
        as.matrix()
}
