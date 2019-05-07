#!/usr/bin/env bash
set -Eeuxo pipefail

# Package checks
#
# See also:
# - R CMD build --help
# - R CMD check --help
# - Travis CI recipe
#   https://github.com/travis-ci/travis-build/blob/master/lib/travis/build/script/r.rb

# Bug fix for `Sys.timezone()` when `timedatectl` is installed.
# https://github.com/rocker-org/rocker-versioned/issues/89
export TZ="America/New_York"

# Get the package version and define the `R CMD build` tarball output.
PKG_NAME="$(basename "$PWD")"
PKG_VERSION="$(grep -E "^Version:\s[.0-9a-z]+$" DESCRIPTION | sed "s/^Version:\s//")"
PKG_TARBALL="${PKG_NAME}_${PKG_VERSION}.tar.gz"

echo "Session information"
Rscript -e "utils::sessionInfo()"
Rscript -e "sessioninfo::session_info()"

echo "R CMD check"
# Set `--as-cran` flag for extra verbose incoming package checks.
R CMD build . --no-build-vignettes --no-manual
R CMD check "$PKG_TARBALL" --ignore-vignettes --no-manual --timings

echo "BiocCheck"
# Note that running `R CMD BiocCheck` directly requires `script/BiocCheck` to
# be installed at `/usr/lib64/R/bin`. Otherwise, can run directly using Rscript.
# Set `--new-package` flag for extra verbose incoming package checks.
Rscript -e "BiocCheck::BiocCheck( \
    package = \"${PKG_TARBALL}\", \
    \`no-check-R-ver\` = TRUE, \
    \`no-check-bioc-help\` = TRUE, \
    \`no-check-remotes\` = TRUE, \
    \`no-check-version-num\` = TRUE, \
    \`no-check-vignettes\` = TRUE, \
    \`quit-with-status\` = TRUE)"

rm "$PKG_TARBALL"

echo "Coverage"
Rscript -e "covr::package_coverage()"

echo "lintr"
Rscript -e "if (packageVersion(\"base\") >= \"3.6\") \
    lintr::lint_package(path = \".\")"
