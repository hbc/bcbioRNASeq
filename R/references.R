#' Download BibTex references
#'
#' If file isn't present, download latest version from the HBC website
#'
#' @author Michael Steinbaugh
#'
#' @param bibtex_file BibTex library file
#'
#' @export
references <- function(bibtex_file = "bcbioRnaseq.bib") {
    if (!file.exists(bibtex_file)) {
        download.file(file.path("http://bioinformatics.sph.harvard.edu",
                                "bcbioRnaseq",
                                bibtex_file),
                      bibtex_file)
    }
}
