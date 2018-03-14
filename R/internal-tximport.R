#' Import RNA-Seq Counts
#'
#' Import RNA-seq counts using [tximport::tximport()].
#'
#' Normalized counts are loaded as length-scaled transcripts per million.
#' https://goo.gl/h6fm15
#'
#' @author Michael Steinbaugh, Rory Kirchner
#' @keywords internal
#' @noRd
#'
#' @importFrom readr read_tsv
#' @importFrom tximport tximport
#'
#' @inheritParams tximport::tximport
#' @param sampleDirs Sample directories to import.
#' @param type *Optional.* Manually specify the expression caller to use.
#'   If `NULL`, if defaults to our preferred priority.
#'
#' @seealso [tximport::tximport()].
#'
#' @return `list` containing count matrices.
.tximport <- function(
    sampleDirs,
    type = NULL,
    txIn = TRUE,
    txOut = FALSE,
    tx2gene
) {
    assert_all_are_dirs(sampleDirs)
    assertIsAStringOrNULL(type)
    if (is_a_string(type)) {
        assert_is_subset(type, validCallers)
    }
    assertIsTx2gene(tx2gene)
    tx2gene <- as.data.frame(tx2gene)

    # Check for count output format, by using the first sample directory
    subdirs <- list.dirs(sampleDirs[[1]], full.names = TRUE, recursive = FALSE)
    assert_are_intersecting_sets(basename(subdirs), validCallers)

    # Set the default priority if caller isn't specified
    if (!is_a_string(type)) {
        if ("salmon" %in% basename(subdirs)) {
            type <- "salmon"
        } else if ("kallisto" %in% basename(subdirs)) {
            type <- "kallisto"
        } else if ("sailfish" %in% basename(subdirs)) {
            type <- "sailfish"
        }
    }

    # Locate `quant.sf` files for salmon or sailfish output
    if (type %in% c("salmon", "sailfish")) {
        files <- list.files(
            path = file.path(sampleDirs, type),
            pattern = "quant.sf",
            full.names = TRUE,
            recursive = TRUE
        )
    } else if (type == "kallisto") {
        files <- list.files(
            path = file.path(sampleDirs, type),
            pattern = "abundance.h5",
            full.names = TRUE,
            recursive = TRUE
        )
    }
    assert_all_are_existing_files(files)
    names(files) <- names(sampleDirs)

    # Begin loading of selected counts
    inform(paste("Reading", type, "counts using tximport"))

    tximport(
        files = files,
        type = type,
        txIn = txIn,
        txOut = txOut,
        countsFromAbundance = "lengthScaledTPM",
        tx2gene = tx2gene,
        ignoreTxVersion = TRUE,
        importer = read_tsv
    )
}



.regenerateTximportList <- function(object) {
    assert_is_all_of(object, "bcbioRNASeq")
    list(
        "abundance" = assays(object)[["tpm"]],
        "counts" = assays(object)[["raw"]],
        "length" = assays(object)[["length"]],
        "countsFromAbundance" = metadata(object)[["countsFromAbundance"]]
    )
}
