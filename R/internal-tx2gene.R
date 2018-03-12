#' Transcript to Gene Annotations
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @importFrom basejump tx2gene
#' @importFrom dplyr mutate
#' @importFrom magrittr set_rownames
#' @importFrom readr read_csv
#'
#' @param projectDir Project directory.
#' @param genomeBuild Genome build.
#' @param release Ensembl release version.
#'
#' @return `data.frame` with unique rownames.
.tx2gene <- function(projectDir, organism, release = NULL) {
    assert_all_are_dirs(projectDir)
    assert_is_a_string(organism)
    assertIsAnImplicitIntegerOrNULL(release)
    inform("Obtaining transcript-to-gene mappings")
    file <- file.path(projectDir, "tx2gene.csv")
    if (file.exists(file)) {
        # bcbio tx2gene
        data <- read_csv(file, col_names = c("txID", "geneID")) %>%
            .[order(.[["txID"]]), ] %>%
            as.data.frame() %>%
            # Ensure transcript versions are removed
            mutate(
                txID = gsub(
                    x = .data[["txID"]],
                    pattern = "\\.\\d+",
                    replacement = ""
                )
            ) %>%
            set_rownames(.[["txID"]])
    } else {
        # Fall back to querying Ensembl
        warn(paste(
            "tx2gene.csv file missing.",
            "Attempting to generate tx2gene from Ensembl instead."
        ))
        data <- tx2gene(organism, release = release)
    }
    assertIsTx2gene(data)
    data
}
