#' Transcript to Gene Annotations
#'
#' @author Michael Steinbaugh
#' @keywords internal
#'
#' @importFrom basejump tx2gene
#' @importFrom dplyr arrange mutate
#' @importFrom magrittr set_rownames
#' @importFrom readr read_csv
#' @importFrom rlang !! sym
#'
#' @param projectDir Project directory.
#' @param genomeBuild Genome build.
#' @param release Ensembl release version.
#'
#' @return [data.frame], with unique rownames.
#' @noRd
.tx2gene <- function(projectDir, organism, release = NULL) {
    assert_all_are_dirs(projectDir)
    assert_is_a_string(organism)
    assert_is_an_implicit_integer_or_null(release)
    file <- file.path(projectDir, "tx2gene.csv")
    if (file.exists(file)) {
        # bcbio tx2gene
        data <- read_csv(file, col_names = c("enstxp", "ensgene")) %>%
            arrange(!!sym("enstxp")) %>%
            as.data.frame() %>%
            # Ensure transcript versions are removed
            mutate(
                enstxp = gsub(
                    x = .data[["enstxp"]],
                    pattern = "\\.\\d+",
                    replacement = "")
            ) %>%
            set_rownames(.[["enstxp"]])
    } else {
        # Fall back to using annotable tx2gene
        warn(paste(
            "tx2gene.csv file missing.",
            "Attempting to generate tx2gene from Ensembl instead."
        ))
        data <- tx2gene(organism, release = release)
    }
    assert_is_tx2gene(data)
    data
}
