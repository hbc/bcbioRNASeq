#' @rdname import
#' @description Import RNA-Seq counts using \code{tximport}
#'
#' @import readr
#' @import stringr
#' @importFrom stats na.omit
#' @importFrom tximport tximport
#' @importFrom utils write.csv
#'
#' @param type The type of software used to generate the abundances. Follows the
#'   conventions of \code{tximport}.
#' @param samples Specify the names of samples in bcbio final directory to
#'   input. If \code{NULL} (default), all samples will be imported.
#'
#' @return txi \code{tximport} list
#' @export
import_counts <- function(
    run,
    type = "salmon",
    grep = NULL,
    samples = NULL) {
    check_run(run)
    if (!type %in% c("salmon", "sailfish")) {
        stop("unsupported counts input format")
    }

    # Sample name grep pattern matching. Run above `samples` to override.
    if (!is.null(grep)) {
        samples <- stringr::str_subset(
            names(run$sample_dirs), pattern = grep)
        if (!length(samples)) {
            stop("grep string didn't match any samples")
        }
    }

    if (is.null(samples)) {
        samples <- names(run$sample_dirs)
    } else {
        samples <- sort(unique(stats::na.omit(samples)))
    }

    # Support for salmon and sailfish file structure
    sample_files <- file.path(run$sample_dirs[samples],
                              type, "quant", "quant.sf")
    names(sample_files) <- names(run$sample_dirs[samples])
    if (!all(file.exists(sample_files))) {
        stop(paste(type, "count files do not exist"))
    }

    # Begin loading of selected counts
    message(paste("loading", type, "counts"))
    message(paste(names(sample_files), collapse = "\n"))
    # writeLines(sample_files)

    tx2gene <- import_file(run, file = "tx2gene.csv", col_names = FALSE)

    # Import the counts
    # https://goo.gl/h6fm15
    # countsFromAbundance = c("no", "scaledTPM", "lengthScaledTPM")
    # Use `lengthScaledTPM` for `salmon` and `sailfish`
    txi <- tximport(files = sample_files,
                    type = type,
                    tx2gene = tx2gene,
                    reader = read_tsv,
                    countsFromAbundance = "lengthScaledTPM")

    return(txi)
}
