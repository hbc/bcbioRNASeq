#!/usr/bin/env Rscript

# Check package lints with lintr.
# Updated 2019-07-15.

options(
    error = quote(quit(status = 1L)),
    warning = quote(quit(status = 1L))
)

requireNamespace("lintr", quietly = TRUE)

lints <- lintr::lint_package(path = ".")

if (length(lints) > 0L) {
    print(lints)
    stop("Package failed lintr checks.")
}
