#' Query the DAVID website
#'
#' @rdname david
#' @author Michael Steinbaugh



#' @rdname david
#' @description Wrapper function that runs gene set enrichment analysis (GSEA)
#'   with the \code{RDAVIDWebService} package, using simplified input options.
#'
#' @param res_tbl \code{res_tables()} list object
#' @param direction Gene direction: up, down, or both
#' @param alpha Alpha level cutoff
#' @param count Minimum hit count
#' @param write_results Write results as tab separated values (TSV) files to
#'   disk
#'
#' @return List of \code{RDAVIDWebService} report objects
#' @export
run_david <- function(
    res_tbl,
    direction = "both",
    alpha = 0.05,
    count = 3,
    write_results = TRUE) {
    # Email
    if (is.null(getOption("email"))) {
        stop("No email found in options().
             Can be saved globally in .Rprofile.")
    }

    # Enrichment direction
    if (direction == "both") {
        foreground <- res_tbl$deg_lfc$ensembl_gene_id
    } else if (direction == "up") {
        foreground <- res_tbl$deg_lfc_up$ensembl_gene_id
    } else if (direction == "down") {
        foreground <- res_tbl$deg_lfc_down$ensembl_gene_id
    } else {
        stop("enrichment direction is required")
    }

    # Count threshold
    if (!is.numeric(count) | length(count) != 1 | count < 0) {
        stop("Please specify the minimum count cutoff of gene hits per
             annotation as a single non-negative numeric")
    }

    name <- res_tbl$name
    contrast <- res_tbl$contrast
    contrast_name <- contrast %>% str_replace_all(" ", "_")

    background <- res_tbl$all$ensembl_gene_id
    id_type = "ENSEMBL_GENE_ID"

    david <- DAVIDWebService$new(
        email = getOption("email"),
        url = "https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")

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
    cutoff_chart <- david %>%
        # S4 DAVIDFunctionalAnnotationChart
        getFunctionalAnnotationChart %>%
        asS3 %>%
        as_tibble %>%
        set_names_snake
    if (nrow(cutoff_chart) > 0) {
        # Filter by gene counts and alpha level
        cutoff_chart <- cutoff_chart %>%
            .[.$count >= count, ] %>%
            .[.$fdr < alpha, ]
        # Set NULL if all processes got filtered
        if (nrow(cutoff_chart) > 0) {
            cutoff_chart <- cutoff_chart %>%
                # FDR should be in 0-1 range, not a percentage!
                # https://goo.gl/Dmf4Qu
                # https://goo.gl/yRPZof
                mutate(fdr = .data$fdr / 100) %>%
                arrange(!!!syms(c("category", "fdr")))
        } else {
            cutoff_chart <- NULL
        }
    } else {
        cutoff_chart <- NULL
    }

    # Write results as TSV files to disk
    if (isTRUE(write_results)) {
        write_dir <- "results/david"
        dir.create(write_dir, recursive = TRUE, showWarnings = FALSE)

        getClusterReportFile(
            david,
            fileName = file.path(
                write_dir,
                paste(name,
                      contrast_name,
                      direction,
                      "cluster_report.tsv",
                      sep = "_")))
        getFunctionalAnnotationChartFile(
            david,
            fileName = file.path(
                write_dir,
                paste(name,
                      contrast_name,
                      direction,
                      "functional_annotation_chart.tsv",
                      sep = "_")))
        getFunctionalAnnotationTableFile(
            david,
            fileName = file.path(
                write_dir,
                paste(name,
                      contrast_name,
                      direction,
                      "functional_annotation_table.tsv",
                      sep = "_")))
        getGeneListReportFile(
            david,
            fileName = file.path(
                write_dir,
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



#' @rdname david
#' @description Generate a DAVID results table for knit report
#'
#' @param david DAVID results list generated by \code{run_david}
#'
#' @return Cutoff chart kable
#' @export
david_table <- function(david) {
    name <- deparse(substitute(david))
    df <- david$cutoff_chart
    if (is.null(df)) {
        stop("DAVID analysis found no significant enrichment")
    }

    df <- df %>%
        as_tibble %>%
        remove_rownames %>%
        .[, c("category", "term", "count", "benjamini")] %>%
        mutate(
            # Add spacing to colon and tilde characters so knitr doesn't output
            # them as links in the report table.
            term = str_replace_all(.data$term, "(~|\\:)", " \\1 "),
            # Truncate the terms so the table isn't too wide
            term = str_trunc(.data$term, side = "right", width = 60),
            # Display BH P values in scientific notation
            # https://goo.gl/Dmf4Qu
            benjamini = format(.data$benjamini, digits = 3, scientific = TRUE))

    # Set the caption
    caption <- paste(name, "david", "functional annotation", sep = " : ")

    return(kable(df, caption = caption))
}
