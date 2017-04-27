#' Utility functions imported from basejump package
#'
#' snake_case variants for function naming consistency
#'
#' @author Michael Steinbaugh



#' @rdname basejump
#' @description Handle multiple kables in a single RMarkdown chunk.
#' Import of \href{https://github.com/steinbaugh/basejump/blob/master/R/kables.R}{kables}.
#'
#' @param list List of tables (e.g. data frame, matrix)
#' @param captions Caption character vector
#'
#' @export
kables <- function(list, captions = NULL) {
    output <- opts_knit$get("rmarkdown.pandoc.to")
    if (!is.null(output)) {
        tables <- lapply(seq_along(list), function(a) {
            kable(list[a], caption = captions[a])
        })
        return(asis_output(tables))
    } else {
        return(list)
    }
}



#' @rdname basejump
#' @description Make syntactically valid names out of character vectors.
#' Variant of \href{https://github.com/steinbaugh/basejump/blob/master/R/makeNames.R}{makeNames}.
#'
#' @param character Character vector
#'
#' @export
make_names <- function(character) {
    character %>%
        as.character %>%
        # Convert non-alphanumeric characters
        gsub("[^[:alnum:]]", ".", .) %>%
        # Combine multiple underscores
        gsub("[\\.]+", ".", .) %>%
        # Strip leading or trailing underscores
        gsub("(^\\.|\\.$)", "", .) %>%
        # Special names
        gsub("(m|nc|r)RNA", "\\1rna", .) %>%
        # Convert acronyms to mixed case
        gsub("([A-Z])([A-Z]+)", "\\1\\L\\2", ., perl = TRUE) %>%
        # Make first letter lowercase
        gsub("(^[A-Z]{1})", "\\L\\1", ., perl = TRUE) %>%
        # Convert camelCase
        gsub("([a-z0-9])([A-Z])", "\\1.\\L\\2", ., perl = TRUE) %>%
        # Ensure syntactically valid names
        make.names %>%
        # Lowercase everything
        tolower
}



#' @rdname basejump
#' @description Make names in snake_case.
#' Variant of \href{https://github.com/steinbaugh/basejump/blob/master/R/makeNames.R}{makeNamesSnake}.
#' @export
make_names_snake <- function(character) {
    character %>% make_names %>% gsub("\\.", "_", .)
}



#' @rdname basejump
#' @description Set names in snake_case.
#' Variant of \href{https://github.com/steinbaugh/basejump/blob/master/R/setNames.R}{setNamesSnake}.
#'
#' @param data Data type that supports name assignments
#'
#' @export
set_names_snake <- function(data) {
    set_names(data, make_names_snake(names(data)))
}
