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
#' @importFrom fs path
#' @importFrom readr read_tsv
#' @importFrom tximport tximport
#'
#' @param sampleDirs Sample directories to import.
#' @param tx2gene Transcript to gene annotations.
#'
#' @seealso
#' - [tximport::tximport()].
#'
#' @return Counts saved in [tximport] list object.
.tximport <- function(sampleDirs, tx2gene) {
    assert_all_are_dirs(sampleDirs)
    assertIsTx2gene(tx2gene)
    tx2gene <- as.data.frame(tx2gene)

    # Check for count output format, by using the first sample directory
    subdirs <- list.dirs(
        sampleDirs[[1]],
        full.names = FALSE,
        recursive = FALSE)

    assert_are_intersecting_sets(c("salmon", "sailfish"), subdirs)
    # Salmon takes priority over sailfish
    if ("salmon" %in% subdirs) {
        type <- "salmon"
    } else if ("sailfish" %in% subdirs) {
        type <- "sailfish"
    }

    # Locate `quant.sf` file for salmon or sailfish output
    if (type %in% c("salmon", "sailfish")) {
        sampleFiles <- list.files(
            path(sampleDirs, type),
            pattern = "quant.sf",
            full.names = TRUE,
            recursive = TRUE)
    }
    assert_all_are_existing_files(sampleFiles)

    # Assign names to sample files
    names(sampleFiles) <- names(sampleDirs)

    # Begin loading of selected counts
    inform(paste("Reading", type, "counts using tximport"))

    # Import the counts (https://goo.gl/h6fm15)
    if (type %in% c("salmon", "sailfish")) {
        countsFromAbundance <- "lengthScaledTPM"
    } else {
        countsFromAbundance <- "no"
    }

    tximport(
        files = sampleFiles,
        type = type,
        tx2gene = tx2gene,
        importer = read_tsv,
        countsFromAbundance = countsFromAbundance,
        ignoreTxVersion = TRUE)
}
