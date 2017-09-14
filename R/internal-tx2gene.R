#' Transcript to Gene Annotations
#'
#' @author Michael Steinbaugh
#' @keywords internal
#'
#' @param projectDir Project directory.
#' @param genomeBuild Genome build.
#' @param release Ensembl release version.
#'
#' @return [data.frame] with unique rownames.
.tx2gene <- function(projectDir, genomeBuild = NULL, release = "current") {
    filePath <- file.path(projectDir, "tx2gene.csv")
    if (file.exists(filePath)) {
        # bcbio tx2gene
        read_csv(filePath, col_names = c("enstxp", "ensgene")) %>%
            arrange(!!sym("enstxp")) %>%
            as.data.frame %>%
            set_rownames(.[["enstxp"]])
    } else {
        # annotable tx2gene
        if (is.null(genomeBuild)) {
            stop("Genome build required for annotable tx2gene")
        }
        tx2gene(genomeBuild, release = release)
    }
}
