#!/usr/bin/env Rscript

## Check package coverage with covr.
## Updated 2019-07-22.

options(
    error = quote(quit(status = 1L)),
    warning = quote(quit(status = 1L))
)

if (!dir.exists("tests")) {
    message("No unit tests defined in `tests/` directory.")
    quit()
}

requireNamespace("covr", quietly = TRUE)

cov <- covr::package_coverage()
pct <- covr::percent_coverage(cov)

if (pct < 90L) {
    print(cov)
    stop(sprintf("Coverage is %s%%.", round(pct, digits = 2L)))
}
