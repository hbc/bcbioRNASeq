#' Transcript to Gene Annotations
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @importFrom basejump tx2gene
#'
#' @param projectDir Project directory.
#' @param organism Organism name.
#' @param genomeBuild Genome build.
#' @param release Ensembl release version.
#'
#' @return `data.frame` with unique rownames.
readTx2gene <- function(file) {
    assert_all_are_dirs(projectDir)
    assertIsAStringOrNULL(organism)
    assertIsAnImplicitIntegerOrNULL(release)
    assertIsAStringOrNULL(genomeBuild)

    inform("Obtaining transcript-to-gene mappings")
    file <- file.path(projectDir, "tx2gene.csv")
    # Warn instead of abort on missing `tx2gene.csv`.
    # Providing fallback support by querying AnnotationHub.
    assert_all_are_existing_files(file, severity = "warning")

    if (file.exists(file)) {
        # bcbio `tx2gene.csv` (recommended)
        data <- read_csv(file, col_names = c("txID", "geneID"))
        # Don't attempt to strip transcript versions
        data <- as.data.frame(data)
        assert_has_no_duplicates(data[["txID"]])
        rownames(data) <- data[["txID"]]
    } else {
        message("Falling back to obtaining mappings from Ensembl")
        assert_is_a_string(organism)
        data <- tx2gene(
            organism = organism,
            release = release,
            genomeBuild = genomeBuild
        )
    }

    assertIsTx2gene(data)
    data
}
