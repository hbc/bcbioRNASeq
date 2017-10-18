.sanitize <- function(tx2gene){
    tx2gene[["enstxp"]] = gsub("\\.[0-9]+", "", tx2gene[["enstxp"]])
    tx2gene
}

#' Transcript to Gene Annotations
#'
#' @author Michael Steinbaugh
#' @keywords internal
#'
#' @param projectDir Project directory.
#' @param genomeBuild Genome build.
#' @param release Ensembl release version.
#'
#' @return [data.frame], with unique rownames.
#' @noRd
.tx2gene <- function(projectDir, organism, release = "current") {
    filePath <- file.path(projectDir, "tx2gene.csv")
    if (file.exists(filePath)) {
        # bcbio tx2gene
        read_csv(filePath, col_names = c("enstxp", "ensgene")) %>%
            arrange(!!sym("enstxp")) %>%
            as.data.frame() %>%
            set_rownames(.[["enstxp"]]) %>%
            .sanitize()
    } else {
        # Fall back to using annotable tx2gene
        warning(paste(
            "tx2gene.csv file missing.",
            "Attempting to generate tx2gene from Ensembl instead."),
            call. = FALSE)
        tx2gene(organism, release = release)
    }
}
