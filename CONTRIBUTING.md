# Contributing to bcbioRNASeq development

- For all changes, fork or create a new branch, then issue a pull request that will be reviewed.
- Do not commit changes directly to master branch.
- Open issues for any changes affecting `bcbioRNASeq` class.
- Support is only provided for the current release version.


## Package style syntax

Attempt to follow these style guides, by priority: (1) [Bioconductor coding style](https://www.bioconductor.org/developers/how-to/coding-style/); (2) [tidyverse style guide](http://style.tidyverse.org); (3) [Google R style guide](https://google.github.io/styleguide/Rguide.xml).

- Object and function parameter arguments should be formated in lowerCamelCase.
- Use spaces instead of tabs. Indent with 4 spaces.
- Use `<-` and not `=` for assignment.
- Explicitly define `if`/`then` statements using opening and closing brackets (`{`, `}`).
- `else` declarations should be on the same line as the closing bracket.
- Use spaces around assignment operators (`<-`), brackets (`{`), and comma delimiters (e.g. `c("foo", "bar")` not `c("foo","bar")`.


## Required checks

If there are changes to `bcbioRNASeq` class or `loadRNASeq()` function, be sure to rebuild all working example data:

```r
source(file.path("data-raw", "examples.R"))
```

For all package updates, run these commands prior to a pull request:

```r
lintr::lint_package()
devtools::document()
devtools::test()
devtools::check()
BiocCheck::BiocCheck(getwd())
```

### Rebuild website

```
unlink("docs", recursive = TRUE)
pkgdown::build_site()
```

Beside this, the package should be able to run the templates. This should be
done in the parent directory of `bcbioRNASeq` folder.

```r
git clone https://github.com/bcbio/bcbio_rnaseq_output_example.git
cd bcbio_rnaseq_output_example
Rscript -e 'testthat::test_file("test_reports.R")'
```

Make sure `qc.html` and `de.html` are created inside `reports` folder and
the figures and tables look as expected.
