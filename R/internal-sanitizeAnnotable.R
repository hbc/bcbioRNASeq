.sanitizeAnnotable <- function(annotable) {
    # Drop any nested list columns (e.g. `entrez`). These's don't play
    # nicely with downstream R Markdown functions.
    nestedCols <- vapply(
        X = annotable,
        FUN = is.list,
        FUN.VALUE = logical(1))
    if (any(nestedCols)) {
        annotable <- annotable[, which(!nestedCols)]
    }
    annotable
}
