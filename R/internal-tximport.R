#' Import RNA-Seq Counts
#'
#' Import RNA-seq counts using [tximport()]. Currently supports
#' [salmon](https://combine-lab.github.io/salmon/) (**recommended**) and
#' [sailfish](http://www.cs.cmu.edu/~ckingsf/software/sailfish/).
#'
#' @author Michael Steinbaugh, Rory Kirchner
#' @keywords internal
#' @noRd
#'
#' @param sampleDirs Sample directories to import.
#' @param tx2gene Transcript to gene annotations.
#'
#' @seealso
#' - [tximport::tximport()].
#'
#' @return Counts saved in [tximport] list object.
.tximport <- function(sampleDirs, tx2gene) {
    # Check for count output format, by using the first sample directory
    subdirs <- list.dirs(sampleDirs[[1]],
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
        sampleFiles <- list.files(
            file.path(sampleDirs, type),
            pattern = "quant.sf",
            full.names = TRUE,
            recursive = TRUE)
    }

    # Assign names to sample files
    names(sampleFiles) <- names(sampleDirs)

    # Begin loading of selected counts
    message(paste("Reading", type, "counts using tximport"))

    # Import the counts (https://goo.gl/h6fm15)
    if (type %in% c("salmon", "sailfish")) {
        countsFromAbundance <- "lengthScaledTPM"
    } else {
        countsFromAbundance <- "no"
    }

    tximport(files = sampleFiles,
             type = type,
             tx2gene = as.data.frame(tx2gene),
             importer = read_tsv,
             countsFromAbundance = countsFromAbundance,
            ignoreTxVersion = TRUE)
}
