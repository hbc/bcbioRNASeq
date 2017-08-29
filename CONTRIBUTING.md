# Contributing to bcbioRNASeq development

- For all changes, fork or create a new branch, then issue a pull request that will be reviewed.
- Do not commit changes directly to master branch.
- Open issues for any changes affecting `bcbioRNADataSet` class.
- Support is only provided for the current release version.


## Required checks

If there are changes to `bcbioRNADataSet` class or `loadRNASeqRun()` function, be sure to rebuild all working example data:

```r
source(file.path("data-raw", "examples.R"))
```

For all package updates, run these commands prior to a pull request:

```r
lintr::lint_package()
devtools::document()
devtools::check()
BiocCheck::BiocCheck(getwd())
pkgdown::build_site()
```
