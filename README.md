# bcbioRnaseq

[![Build Status](https://travis-ci.org/hbc/bcbioRnaseq.svg?branch=master)](https://travis-ci.org/hbc/bcbioRnaseq)
[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![codecov](https://codecov.io/gh/hbc/bcbioRnaseq/branch/master/graph/badge.svg)](https://codecov.io/gh/hbc/bcbioRnaseq)

Quality control and differential expression for [bcbio][] RNA-seq experiments.


## Installation

This is an [R][] package.

### [Bioconductor][] method

```r
source("https://bioconductor.org/biocLite.R")
biocLite("hbc/bcbioRnaseq")
```

### [devtools][] method

```r
install.packages("devtools")
devtools::install_github("hbc/bcbioRnaseq")
```


## Examples

[HTML reports](http://bcb.io/bcbio_rnaseq_output_example) rendered from the default RMarkdown templates included in the package.

- [Quality control](http://bcb.io/bcbio_rnaseq_output_example/qc.html)
- [Differential expression](http://bcb.io/bcbio_rnaseq_output_example/de.html)


[bcbio]: https://github.com/chapmanb/bcbio-nextgen
[bioconductor]: https://bioconductor.org
[devtools]: https://cran.r-project.org/package=devtools
[r]: https://www.r-project.org
[rmarkdown]: http://rmarkdown.rstudio.com
