#' Print top table
#'
#' @author Michael Steinbaugh
#'
#' @importFrom utils head
#'
#' @param df Data frame
#' @param caption Caption
#'
#' @return Data frame or kable, depending on the call
#' @export
top_table <- function(df, caption = NULL) {
    discard <- c("description",
                 "lfc_se",
                 "pvalue",
                 "stat")
    df %>%
        .[, setdiff(colnames(.), discard)] %>%
        utils::head(n = 50) %>%
        print_table(caption = caption)
}
