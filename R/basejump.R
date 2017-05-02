#' Utility functions from basejump package
#'
#' @rdname basejump
#'
#' @author Michael Steinbaugh
#' 
#' @param dna DNA sequence (ATGC nucleotides)
#' @param list List of tables (e.g. data frame, matrix)
#' @param captions Caption character vector
#' @param character Character vector
#' @param dir Output directory
#' @param data Data type that supports name assignments



## https://github.com/steinbaugh/basejump/blob/dev/R/dna.R
#' @rdname basejump
#' @keywords internal
#' @description Complement DNA sequence
#' @export
comp <- function(dna) {
    dna <- toupper(dna)
    comp <- dna %>%
        # AT base pair swap
        gsub("A", "A1", .) %>%
        gsub("T", "A", .) %>%
        gsub("A1", "T", .) %>%
        # GC base pair swap
        gsub("G", "G1", .) %>%
        gsub("C", "G", .) %>%
        gsub("G1", "C", .)
    return(comp)
}



## https://github.com/steinbaugh/basejump/blob/dev/R/kables.R
#' @rdname basejump
#' @keywords internal
#' @description Handle multiple kables in a single RMarkdown chunk
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



## https://github.com/steinbaugh/basejump/blob/dev/R/makeNames.R
#' @rdname basejump
#' @keywords internal
#' @description Make syntactically valid names out of character vectors
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



## https://github.com/steinbaugh/basejump/blob/dev/R/makeNames.R
#' @rdname basejump
#' @keywords internal
#' @description Make names in snake_case
#' @export
make_names_snake <- function(character) {
    character %>% make_names %>% gsub("\\.", "_", .)
}



## https://github.com/steinbaugh/basejump/blob/dev/R/dna.R
#' @rdname basejump
#' @keywords internal
#' @description Reverse complement DNA sequence
#' @export
revcomp <- function(dna) {
    dna <- toupper(dna)
    comp <- comp(dna)
    revcomp <- strsplit(comp, "")[[1]] %>%
        .[order(seq_along(.), decreasing = TRUE)] %>%
        paste(., sep = "", collapse = "")
    return(revcomp)
}



## https://github.com/steinbaugh/basejump/blob/dev/R/saveData.R
#' @rdname basejump
#' @keywords internal
#' @description Quick save to \code{data} directory
#' @export
save_data <- function(..., dir = "data") {
    if (!isString(dir)) {
        stop("dir must be a string")
    }
    if (!dir.exists(dir)) {
        dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    }
    names <- as.character(substitute(list(...)))[-1L]
    objs <- get_objs_from_dots(dots(...))
    paths <- file.path(dir, paste0(objs, ".rda"))
    message(paste("Saving", toString(names), "to", dir))
    mapply(save, list = objs, file = paths)
    invisible()
}



#' @rdname basejump
#' @description Set names in dot notation
#' @export
set_names_dot <- function(data) {
    set_names(data, make_names(colnames(data)))
}



## https://github.com/steinbaugh/basejump/blob/dev/R/setNames.R
#' @rdname basejump
#' @keywords internal
#' @description Set names in snake_case
#' @export
set_names_snake <- function(data) {
    set_names(data, make_names_snake(names(data)))
}






# Unexported functions ====
# https://github.com/hadley/devtools/blob/dev/R/utils.r
dots <- function(...) {
    eval(substitute(alist(...)))
}



# https://github.com/hadley/devtools/blob/dev/R/infrastructure.R
get_objs_from_dots <- function(.dots) {
    if (length(.dots) == 0L) {
        stop("Nothing to save", call. = FALSE)
    }
    is_name <- vapply(.dots, is.symbol, logical(1))
    if (any(!is_name)) {
        stop("Can only save existing named objects", call. = FALSE)
    }
    objs <- vapply(.dots, as.character, character(1))
    duplicated_objs <- which(set_names(duplicated(objs), objs))
    if (length(duplicated_objs) > 0L) {
        objs <- unique(objs)
        warning("Saving duplicates only once: ",
                paste(names(duplicated_objs), collapse = ", "),
                call. = FALSE)
    }
    return(objs)
}
