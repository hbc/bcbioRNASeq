#' Prepare `SummarizedExperiment()`
#'
#' This is a utility wrapper for `SummarizedExperiment()` that provides
#' automatic subsetting for `colData` and `rowData`.
#'
#' This function also provides automatic metadata slotting of multiple useful
#' environment parameters:
#'
#' - `date`: Today's date.
#' - `wd`: Working directory.
#' - `utilsSessionInfo`: [utils::sessionInfo()] return.
#' - `devtoolsSessionInfo`: [devtools::session_info()] return.
#'
#' @family bcbio Utilities
#' @keywords internal
#'
#' @importFrom devtools session_info
#' @importFrom magrittr set_rownames
#' @importFrom S4Vectors head
#' @importFrom scales percent
#' @importFrom tibble column_to_rownames has_rownames
#' @importFrom utils sessionInfo
#'
#' @param assays List containing RNA-seq count matrices with matching
#'   dimensions. Counts can be passed in either dense (`matrix`) or sparse
#'   (`dgCMatrix`, `dgTMatrix`) format.
#' @param rowData Object describing assay matrix rows. Must support
#'   [base::dim()].
#' @param colData Object describing assay matrix columns. Must support
#'   [base::dim()].
#' @param metadata *Optional*. Metadata list.
#'
#' @seealso
#' - [SummarizedExperiment::SummarizedExperiment()].
#' - [base::Sys.Date()].
#' - [base::getwd()].
#' - [utils::sessionInfo()].
#'
#' @return [SummarizedExperiment].
#' @export
#'
#' @examples
#' mat <- matrix(
#'     seq(1L:16L),
#'     nrow = 4L,
#'     ncol = 4L,
#'     dimnames = list(
#'         c("ENSMUSG00000000001",
#'           "ENSMUSG00000000003",
#'           "ENSMUSG00000000028",
#'           "ENSMUSG00000000031"),
#'         c("sample_1",
#'           "sample_2",
#'           "sample_3",
#'           "sample_4")))
#' assays <- list(assay = mat)
#' rowData <- data.frame(
#'     ensgene = c("ENSMUSG00000000001",
#'                 "ENSMUSG00000000003",
#'                 "ENSMUSG00000000028",
#'                 "ENSMUSG00000000031"),
#'     biotype = c("coding",
#'                 "coding",
#'                 "coding",
#'                 "pseudogene"),
#'     row.names = rownames(mat))
#' colData <- data.frame(
#'     genotype = c("wildtype",
#'                  "wildtype",
#'                  "knockout",
#'                  "knockout"),
#'     age = c(3L, 6L, 3L, 6L),
#'     row.names = colnames(mat))
#' .prepareSummarizedExperiment(
#'     assays = assays,
#'     rowData = rowData,
#'     colData = colData)
.prepareSummarizedExperiment <- function(
    assays,
    rowData,
    colData,
    metadata = NULL) {
    # Assays ====
    if (!is.list(assays)) {
        stop("'assays' must be a list", call. = FALSE)
    }
    assay <- assays[[1L]]
    if (is.null(dim(assay))) {
        stop("Assay object must support 'dim()'", call. = FALSE)
    }
    # Check for potential dimnames problems
    if (is.null(rownames(assay))) {
        stop("Assay missing rownames", call. = FALSE)
    }
    if (is.null(colnames(assay))) {
        stop("Assay missing colnames", call. = FALSE)
    }
    if (any(duplicated(rownames(assay)))) {
        stop("Non-unique rownames", call. = FALSE)
    }
    if (any(duplicated(colnames(assay)))) {
        stop("Non-unique colnames", call. = FALSE)
    }
    if (!identical(
        make.names(rownames(assay), unique = TRUE, allow_ = TRUE),
        rownames(assay)
    )) {
        stop(paste(
            "Rownames are not valid.",
            "See 'base::make.names()' for more information."
        ), call. = FALSE)
    }
    if (!identical(
        make.names(colnames(assay), unique = TRUE, allow_ = TRUE),
        colnames(assay)
    )) {
        stop(paste(
            "Colnames are not valid.",
            "See 'base::make.names()' for more information."
        ), call. = FALSE)
    }

    # rowData ====
    if (is.null(dim(rowData))) {
        stop("rowData must support 'dim()'", call. = FALSE)
    }
    rowData <- as.data.frame(rowData)
    # Handle tibble rownames
    if (!has_rownames(rowData) &
        "rowname" %in% colnames(rowData)) {
        rowData <- column_to_rownames(rowData)
    }
    if (!has_rownames(rowData)) {
        stop("rowData missing rownames", call. = FALSE)
    }
    if (!all(rownames(assay) %in% rownames(rowData))) {
        missingGenes <- setdiff(rownames(assay), rownames(rowData))
        # Warn instead of stop here, for more lenient handling of deprecated
        # identifiers in a newer Ensembl release
        warning(paste(
            "Ensembl retired genes detected",
            paste0("(", pct(length(missingGenes) / nrow(assay)), ")")
        ), call. = FALSE)
    } else {
        missingGenes <- NULL
    }

    # Fit the annotable annotations to match the number of rows
    rowData <- rowData[rownames(assay), , drop = FALSE]
    rownames(rowData) <- rownames(assay)
    # Ensure all identifiers are set in ensgene, even those without current
    # Ensembl annotation information
    rowData[["ensgene"]] <- rownames(assay)
    rowData <- as(rowData, "DataFrame")

    # colData ====
    if (is.null(dim(colData))) {
        stop("colData must support 'dim()'", call. = FALSE)
    }
    colData <- as.data.frame(colData)
    # Handle tibble rownames
    if (!has_rownames(colData) &
        "rowname" %in% colnames(colData)) {
        colData <- column_to_rownames(colData)
    }
    if (!has_rownames(colData)) {
        stop("colData missing rownames", call. = FALSE)
    }
    if (!all(colnames(assay) %in% rownames(colData))) {
        missingSamples <- setdiff(colnames(assay), rownames(colData))
        stop(paste(
            "Sample mismatch between 'assay' and 'colData':",
            toString(head(missingSamples))
        ), call. = FALSE)
    }
    colData <- colData %>%
        .[colnames(assay), , drop = FALSE] %>%
        set_rownames(colnames(assay)) %>%
        as("DataFrame")

    # Metadata ====
    if (is.null(metadata)) {
        metadata <- list()
    } else {
        if (!is(metadata, "list")) {
            stop("Metadata must be 'list' class object", call. = FALSE)
        }
    }
    metadata[["date"]] <- Sys.Date()
    metadata[["wd"]] <- getwd()
    metadata[["utilsSessionInfo"]] <- utils::sessionInfo()
    metadata[["devtoolsSessionInfo"]] <- devtools::session_info()
    metadata[["missingGenes"]] <- missingGenes

    SummarizedExperiment(
        assays = assays,
        rowData = rowData,
        colData = colData,
        metadata = metadata)
}
