#' Run DAVID analysis
#'
#' Wrapper function that performs gene set enrichment analysis (GSEA) with the
#' \code{RDAVIDWebService} package, using simplified input options.
#'
#' @author Michael Steinbaugh
#'
#' @import dplyr
#' @import readr
#' @import stringr
#' @importFrom RDAVIDWebService addList DAVIDWebService getAnnotationSummary
#'   getClusterReport getClusterReportFile getFunctionalAnnotationChart
#'   getFunctionalAnnotationChartFile getFunctionalAnnotationTable
#'   getFunctionalAnnotationTableFile getGeneCategoriesReport getGeneListReport
#'   getGeneListReportFile setTimeOut
#'
#' @param res_tables \code{res_tables()} list object
#' @param direction Gene direction: up, down, or both
#' @param write_files Write files to disk (\code{TRUE/FALSE})
#' @param count Minimum hit count
#'
#' @return List of \code{RDAVIDWebService} report objects
#' @export
run_david <- function(
    res_tables,
    direction = "both",
    count = 3,
    write_files = TRUE) {
    # Email
    if (is.null(getOption("email"))) {
        stop("no email found in options().
             can be saved globally in .Rprofile.")
    }

    # Enrichment direction
    if (direction == "both") {
        foreground <- res_tables$deg_lfc$ensembl_gene_id
    } else if (direction == "up") {
        foreground <- res_tables$deg_lfc_up$ensembl_gene_id
    } else if (direction == "down") {
        foreground <- res_tables$deg_lfc_down$ensembl_gene_id
    } else {
        stop("enrichment direction is required")
    }

    # Count threshold
    if (!is.numeric(count) | length(count) != 1 | count < 0) {
        stop("please specify the minimum count cutoff of gene hits per
             annotation as a single non-negative numeric")
    }

    name <- res_tables$name
    contrast <- res_tables$contrast
    contrast_name <- contrast %>% gsub(" ", "_", .)

    background <- res_tables$all$ensembl_gene_id
    fdr <- res_tables$alpha
    id_type = "ENSEMBL_GENE_ID"

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
                  "fdr")] %>%
            # FDR should be presented on 0-1 scale, not as a percentage
            mutate_(.dots = set_names(list(~fdr / 100), "fdr")) %>%
            arrange_(.dots = c("category", "fdr")) %>%
            .[.$count >= count, ] %>%
            .[.$fdr < fdr, ]
        # Set NULL if everything got filtered
        if (nrow(cutoff_chart) == 0) {
            cutoff_chart <- NULL
        }
    } else {
        cutoff_chart <- NULL
    }

    # Save the TSV files to disk
    if (isTRUE(write_files)) {
        write_dir <- "results/david"
        dir.create(write_dir, recursive = TRUE, showWarnings = FALSE)

        getClusterReportFile(
            david,
            fileName = file.path(write_dir,
                                 paste(name,
                                       contrast_name,
                                       direction,
                                       "cluster_report.tsv",
                                       sep = "_")))
        getFunctionalAnnotationChartFile(
            david,
            fileName = file.path(write_dir,
                                 paste(name,
                                       contrast_name,
                                       direction,
                                       "functional_annotation_chart.tsv",
                                       sep = "_")))
        getFunctionalAnnotationTableFile(
            david,
            fileName = file.path(write_dir,
                                 paste(name,
                                       contrast_name,
                                       direction,
                                       "functional_annotation_table.tsv",
                                       sep = "_")))
        getGeneListReportFile(
            david,
            fileName = file.path(write_dir,
                                 paste(name,
                                       contrast_name,
                                       direction,
                                       "gene_list_report_file.tsv",
                                       sep = "_")))

        if (!is.null(cutoff_chart)) {
            write_tsv(
                cutoff_chart,
                file.path(write_dir,
                          paste(name,
                                contrast_name,
                                direction,
                                "cutoff_chart.tsv",
                                sep = "_"))
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
