#' Create project directory structure
#'
#' @author Michael Steinbaugh
#' @keywords internal
#'
#' @export
create_project <- function() {
    dir.create("data", showWarnings = FALSE)
    dir.create("figures", showWarnings = FALSE)
    dir.create("meta", showWarnings = FALSE)
    dir.create("results", showWarnings = FALSE)
}
