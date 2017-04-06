#' Create project directory structure
#'
#' @author Michael Steinbaugh
#' @keywords internal
#'
#' @export
create_local_project <- function() {
    local_dirs <- c("cache",
                    "data",
                    "figures",
                    "meta",
                    "results")
    lapply(seq_along(local_dirs), function(a) {
        dir.create(local_dirs[a], showWarnings = FALSE)
    }) %>% invisible
}
