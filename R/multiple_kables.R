#' Knit multiple kables in a single chunk
#'
#' @author Michael Steinbaugh
#' @keywords internal
#'
#' @importFrom knitr asis_output
#'
#' @param list List of tables (e.g. data frame, matrix)
#'
#' @return Knit tables, using \code{kable()}
#' @export
#'
#' @examples
#' multiple_kables(list(iris, mtcars))
multiple_kables <- function(list) {
    tables <- lapply(seq_along(list), function(a) {
        kable(list[a], caption = names(list)[a])
    })
    return(knitr::asis_output(tables))
}
