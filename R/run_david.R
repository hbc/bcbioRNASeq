#' Run DAVID analysis
#'
#' Wrapper function that performs gene set enrichment analysis (GSEA) with the
#' \code{RDAVIDWebService} package, using simplified input options.
#'
#' @author Michael Steinbaugh
#'
#' @import dplyr
#' @import readr
#' @importFrom RDAVIDWebService addList DAVIDWebService getAnnotationSummary
#'   getClusterReport getClusterReportFile getFunctionalAnnotationChart
#'   getFunctionalAnnotationChartFile getFunctionalAnnotationTable
#'   getFunctionalAnnotationTableFile getGeneCategoriesReport getGeneListReport
#'   getGeneListReportFile setTimeOut
#'
#' @param foreground Foreground identifiers
#' @param background Background identifiers
#' @param id_type Identifier type (see DAVID website)
#' @param write_files Write files to disk (\code{TRUE/FALSE})
#' @param write_prefix Add an optional prefix to the file names
#' @param write_dir Directory to write the TSV files
#' @param count Minimum hit count
#' @param alpha False discovery rate cutoff (alpha)
#'
#' @return List of \code{RDAVIDWebService} report objects
#' @export
run_david <- function(
    foreground,
    background = NULL,
    id_type = "ENSEMBL_GENE_ID",
    write_files = TRUE,
    write_prefix = NULL,
    write_dir = "results/david",
    count = 3,
    alpha = 0.1) {
    if (is.null(getOption("email"))) {
        stop("no email found in options().
             can be saved globally in .Rprofile.")
    }
    if (is.null(foreground)) {
        stop("A foreground gene vector is required.")
    }
    if (is.null(id_type) | length(id_type) != 1) {
        stop("identifier string is required")
    }
    if (!is.numeric(count) | length(count) != 1 | count < 0) {
        stop("please specify the minimum count cutoff of gene hits per
             annotation as a single non-negative numeric")
    }
    if (!is.numeric(alpha) | length(alpha) != 1 | alpha < 0 | alpha > 1) {
        stop("please specify the false discovery rate (FDR) cutoff as a single
             numeric in the range of 0-1")
    }

    david <- DAVIDWebService$new(
        email = getOption("email"),
        url = "https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/"
    )

    # Set a longer timeout (30000 default)
    setTimeOut(david, 200000)
    # getTimeOut(david)

    addList(david,
            foreground,
            idType = id_type,
            listName = "Gene",
            listType = "Gene")

    # Set the background only for screens and microarrays
    # Use the organism default for RNA-Seq
    if (!is.null(background)) {
        addList(david,
                background,
                idType = id_type,
                listName = "Background",
                listType = "Background")
    }

    # Generate the annotation chart with cutoffs applied
    cutoff_chart <- getFunctionalAnnotationChart(david) %>% as.data.frame
    if (nrow(cutoff_chart) > 0) {
        cutoff_chart <- cutoff_chart %>%
            set_names_snake %>%
            .[, c("category",
                  "term",
                  "count",
                  "genes",
                  "pvalue",
                  "alpha")] %>%
            dplyr::rename_(.dots = c("p" = "pvalue")) %>%
            # FDR should be presented on 0-1 scale, not as a percentage
            dplyr::mutate_(.dots = set_names(list(~fdr / 100), "fdr")) %>%
            dplyr::arrange_(.dots = c("category", "fdr")) %>%
            .[.$count >= count, ] %>%
            .[.$fdr < alpha, ]
        # Set NULL if everything got filtered
        if (nrow(cutoff_chart) == 0) {
            cutoff_chart <- NULL
        }
    } else {
        cutoff_chart <- NULL
    }

    # Save the TSV files to disk
    if (isTRUE(write_files)) {
        dir.create(write_dir,
                   recursive = TRUE,
                   showWarnings = FALSE)

        if (!is.null(write_prefix)) {
            prefix <- paste0(write_prefix, "_")
        } else {
            prefix <- ""
        }

        getClusterReportFile(
            david,
            fileName = file.path(write_dir,
                                 paste0(prefix,
                                        "cluster_report.tsv")))
        getFunctionalAnnotationChartFile(
            david,
            fileName = file.path(write_dir,
                                 paste0(prefix,
                                        "functional_annotation_chart.tsv")))
        getFunctionalAnnotationTableFile(
            david,
            fileName = file.path(write_dir,
                                 paste0(prefix,
                                        "functional_annotation_table.tsv")))
        getGeneListReportFile(
            david,
            fileName = file.path(write_dir,
                                 paste0(prefix,
                                        "gene_list_report_file.tsv")))

        if (!is.null(cutoff_chart)) {
            readr::write_tsv(
                cutoff_chart,
                file.path(write_dir,
                          paste0(prefix, "cutoff_chart.tsv"))
            )
        }
    }

    # Package DAVID output into a list
    return(list(
        annotation_summary = getAnnotationSummary(david),
        cluster_report = getClusterReport(david),
        cutoff_chart = cutoff_chart,
        functional_annotation_chart = getFunctionalAnnotationChart(david),
        functional_annotation_table = getFunctionalAnnotationTable(david),
        gene_categories_report = getGeneCategoriesReport(david),
        gene_list_report = getGeneListReport(david)
    ))
}
