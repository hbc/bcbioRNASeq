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
biocLite("hbc/bcbioRNASeq")
```

### [devtools][] method

```r
install.packages("devtools")
devtools::install_github("hbc/bcbioRNASeq")
```


## Load [bcbio][] run output

```r
library(bcbioRNASeq)
bcb <- loadRNASeqRun(
    file.path("bcbio_rnaseq_run", "final"),
    interestingGroups = c("genotype", "treatment"))
```

Parameters:

- Primary object: Directory path to the final [bcbio][] run output.
- `interestingGroups`: Character vector with the variables to use for representation of the data that should match columns of interest in the sample metadata.

Consult `help("loadRNASeqRun", "bcbioRNASeq")` for additional documentation.


## [RMarkdown][] templates

This package provides multiple [RMarkdown][] templates, including Quality Control and Differential Expression using DESeq2, which are available in [RStudio][] at `File -> New File -> R Markdown... -> From Template`.

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
