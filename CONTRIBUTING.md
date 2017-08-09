# Contributing to bcbioRnaseq development

- For all changes, fork or create a new branch, then issue a pull request that will be reviewed.
- Do not commit changes directly to master branch.
- Open issues for any changes affecting `bcbioRNADataSet` class.
- Support is only provided for the current release version.


## Required checks

```r
lintr::lint_package()
devtools::document()
devtools::check()
BiocCheck::BiocCheck(getwd())
pkgdown::build_site()
```
