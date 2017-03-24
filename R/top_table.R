#' Print top table
#'
#' @author Michael Steinbaugh
#'
#' @importFrom basejump printTable
#'
#' @param df Data frame
#' @param caption Caption
#'
#' @return Data frame or kable, depending on the call
#' @export
top_table <- function(df, caption = NULL) {
    df %>%
        dplyr::select(-c(description, lfc_se, pvalue, stat)) %>%
        head(n = 50) %>%
        basejump::printTable(caption = caption)
}
