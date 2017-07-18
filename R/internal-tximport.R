#' Import RNA-seq counts
#'
#' Import RNA-seq counts using [tximport()]. Currently supports
#' [salmon](https://combine-lab.github.io/salmon/) (**recommended**) and
#' [sailfish](http://www.cs.cmu.edu/~ckingsf/software/sailfish/).
#'
#' @rdname tximport
#' @keywords internal
#'
#' @author Michael Steinbaugh
#' @author Rory Kirchner
#'
#' @param sample_dirs Sample directories to import.
#' @param tx2gene Transcript to gene annotations.
#'
#' @return Counts saved in [tximport] list object.
#'
#' @seealso
#' - [tximport::tximport()].
.tximport <- function(sample_dirs, tx2gene) {
    # Check for count output format, by using the first sample directory
    subdirs <- list.dirs(sample_dirs[[1L]],
                         full.names = FALSE,
                         recursive = FALSE)
    if ("salmon" %in% subdirs) {
        type <- "salmon"
    } else if ("sailfish" %in% subdirs) {
        type <- "sailfish"
    } else {
        stop("Unsupported counts output format")
    }

    # Locate `quant.sf` file for salmon or sailfish output
    if (type %in% c("salmon", "sailfish")) {
        sample_files <- list.files(
            file.path(sample_dirs, type),
            pattern = "quant.sf",
            full.names = TRUE,
            recursive = TRUE)
    }

    # Assign names to sample files
    names(sample_files) <- names(sample_dirs)

    # Begin loading of selected counts
    message(paste("Reading", type, "counts"))

    # Import the counts (https://goo.gl/h6fm15)
    if (type %in% c("salmon", "sailfish")) {
        counts_from_abundance <- "lengthScaledTPM"
    } else {
        counts_from_abundance <- "no"
    }

    tximport(files = sample_files,
             type = type,
             tx2gene = as.data.frame(tx2gene),
             importer = read_tsv,
             countsFromAbundance = counts_from_abundance)
}
