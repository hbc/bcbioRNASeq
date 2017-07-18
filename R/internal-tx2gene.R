#' Read transcript to gene (tx2gene) annotation file
#'
#' @rdname tx2gene
#' @keywords internal
#' @author Michael Steinbaugh
#'
#' @param project_dir Project directory.
#' @param genome_build Genome build.
.tx2gene <- function(project_dir, genome_build) {
    file_path <- file.path(project_dir, "tx2gene.csv")
    if (file.exists(file_path)) {
        # bcbio-nextgene tx2gene
        read_csv(file_path, col_names = c("enstxp", "ensgene")) %>%
            arrange(!!sym("enstxp")) %>%
            as.data.frame %>%
            set_rownames(.[["enstxp"]])
    } else {
        # annotable tx2gene
        if (is.null(genome_build)) {
            stop("Genome build required for annotable tx2gene")
        }
        tx2gene(genome_build)
    }
}
