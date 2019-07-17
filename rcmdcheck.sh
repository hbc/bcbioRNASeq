#!/usr/bin/env bash
set -Eeuxo pipefail

# R package checks
# Updated 2019-07-16.
#
# See also:
# - R CMD build --help
# - R CMD check --help
# - Travis CI recipe
#   https://github.com/travis-ci/travis-build/blob/master/lib/travis/build/script/r.rb

# Don't require suggested packages to be installed.
export _R_CHECK_FORCE_SUGGESTS_=false

# Bug fix for `Sys.timezone()` when `timedatectl` is installed.
# https://github.com/rocker-org/rocker-versioned/issues/89
export TZ="America/New_York"

# Get the package version and define the `R CMD build` tarball output.
PKG_NAME="$(basename "$PWD")"
PKG_VERSION="$(grep -E "^Version:\s[.0-9a-z]+$" DESCRIPTION | sed "s/^Version:[[:space:]]//")"
PKG_TARBALL="${PKG_NAME}_${PKG_VERSION}.tar.gz"

echo "travis_fold:start:session_info"
Rscript -e "utils::sessionInfo()"
Rscript -e "sessioninfo::session_info()"
echo "travis_fold:end:session_info"

echo "travis_fold:start:r_cmd_check"
# Set `--as-cran` flag for extra verbose incoming package checks.
R CMD build . --no-build-vignettes --no-manual
R CMD check "$PKG_TARBALL" --ignore-vignettes --no-manual --timings
echo "travis_fold:end:r_cmd_check"

echo "travis_fold:start:bioc_check"
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
echo "travis_fold:end:bioc_check"

rm "$PKG_TARBALL"

echo "travis_fold:start:lints"
./lints.R
echo "travis_fold:end:lints"

echo "travis_fold:start:coverage"
./coverage.R
echo "travis_fold:end:coverage"
