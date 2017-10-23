# bcbioRNASeq

[![Build Status](https://travis-ci.org/hbc/bcbioRNASeq.svg?branch=master)](https://travis-ci.org/hbc/bcbioRNASeq)
[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![codecov](https://codecov.io/gh/hbc/bcbioRNASeq/branch/master/graph/badge.svg)](https://codecov.io/gh/hbc/bcbioRNASeq)

Quality control and differential expression for [bcbio][] RNA-seq experiments.


## Installation

This is an [R][] package.

### [Bioconductor][] method

```r
source("https://bioconductor.org/biocLite.R")
```

Ensure that these dependencies install correctly:

```r
biocLite(c(
    "steinbaugh/basejump",
    "lpantano/DEGreport",
    "tidyverse",
    "rmarkdown",
    "knitr",
    "formatR",
    "GenomeInfoDbData"
))
```

Now you're ready to install bcbioRNASeq:

```r
biocLite("hbc/bcbioRNASeq")
```

For the *Functional Analysis* R Markdown template, these additional libraries are required:

```r
biocLite(c(
    "clusterProfiler",
    "DOSE",
    "pathview"
))
```


## Load [bcbio][] run

```r
library(bcbioRNASeq)
bcb <- loadRNASeq(
    uploadDir = file.path("bcbio_rnaseq_run", "final"),
    interestingGroups = c("genotype", "treatment"))
```

This will return a `bcbioRNASeq` object, which is an extension of the [Bioconductor][] [SummarizedExperiment][] container class.

Parameters:

- `uploadDir`: Path to the [bcbio][] final upload directory.
- `interestingGroups`: Character vector of the column names of interest in the sample metadata, which is stored in the `colData()` accessor slot of the `bcbioRNASeq` object. These values should be formatted in camelCase, and can be reassigned in the object after creation (e.g. `metadata(bcb)$interestingGroups <- c("batch", "age")`). They are used for data visualization in the quality control utility functions.

Consult `help("loadRNASeq", "bcbioRNASeq")` for additional documentation.


## [RMarkdown][] templates

This package provides multiple [RMarkdown][] templates, including Quality Control and Differential Expression using [DESeq2][], which are available in [RStudio][] at `File` -> `New File` -> `R Markdown...` -> `From Template`.

### Examples

View example [HTML reports](http://bcb.io/bcbio_rnaseq_output_example) rendered from the default [RMarkdown][] templates included in the package:

- [Quality Control](http://bcb.io/bcbio_rnaseq_output_example/qc-master.html)
- [Differential Expression](http://bcb.io/bcbio_rnaseq_output_example/de-master.html)


[bcbio]: https://github.com/chapmanb/bcbio-nextgen
[Bioconductor]: https://bioconductor.org
[DESeq2]: https://bioconductor.org/packages/release/bioc/html/DESeq2.html
[devtools]: https://cran.r-project.org/package=devtools
[R]: https://www.r-project.org
[RMarkdown]: http://rmarkdown.rstudio.com
[RStudio]: https://www.rstudio.com
[SummarizedExperiment]: http://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html
