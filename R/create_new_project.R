#' Create new project directory structure
#'
#' @author Michael Steinbaugh
#'
#' @export
create_new_project <- function() {
    local_dirs <- c("data",
                    "docs",
                    "figures",
                    "meta",
                    "results")
    lapply(seq_along(local_dirs), function(a) {
        dir.create(local_dirs[a], showWarnings = FALSE)
    }) %>% invisible
}
