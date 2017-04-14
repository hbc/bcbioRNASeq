[Bioconductor]: https://bioconductor.org
[R]: https://www.r-project.org

[`devtools`]: https://cran.r-project.org/package=devtools
[`bcbio-nextgen`]: https://bcbio-nextgen.readthedocs.io



# bcbioRnaseq

[![Project Status: WIP - Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)

Quality control and differential expression for [`bcbio-nextgen`][] RNA-Seq experiments


## Installation

This is an [R][] package.

### [Bioconductor][] method

```{r}
source("https://bioconductor.org/biocLite.R")
biocLite("roryk/bcbioSinglecell")
```

### [`devtools`][] method

```{r}
install.packages("devtools")
devtools::install_github("roryk/bcbioSinglecell")
```
