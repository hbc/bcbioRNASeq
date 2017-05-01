[Bioconductor]: https://bioconductor.org
[bcbio-nextgen]: https://github.com/chapmanb/bcbio-nextgen
[devtools]: https://cran.r-project.org/package=devtools
[R]: https://www.r-project.org



# bcbioRnaseq

[![Build Status](https://travis-ci.org/roryk/bcbioRnaseq.svg?branch=master)](https://travis-ci.org/roryk/bcbioRnaseq)
[![Project Status: WIP - Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)

Quality control and differential expression for [bcbio-nextgen][] RNA-seq experiments.


## Installation

This is an [R][] package.

### [Bioconductor][] method

```{r}
source("https://bioconductor.org/biocLite.R")
biocLite("roryk/bcbioRnaseq")
```

### [devtools][] method

```{r}
install.packages("devtools")
devtools::install_github("roryk/bcbioRnaseq")
```
