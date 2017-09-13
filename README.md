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

## Load bcbio results

```
library(bcbioRNASeq)
bcb <- loadRNASeqRun(
params$uploadDir,
interestingGroups = params$interestingGroups,
params = params)
saveData(bcb, dir = dataDir)
```

Use the Quality Control or/and the Differential Expression template available
in Rstudio at `File -> New File -> R markdown -> From template`.


## Examples

[HTML reports](http://bcb.io/bcbio_rnaseq_output_example) rendered from the default RMarkdown templates included in the package.

- [Quality control](http://bcb.io/bcbio_rnaseq_output_example/qc-master.html)
- [Differential expression](http://bcb.io/bcbio_rnaseq_output_example/de-master.html)


[bcbio]: https://github.com/chapmanb/bcbio-nextgen
[Bioconductor]: https://bioconductor.org
[devtools]: https://cran.r-project.org/package=devtools
[R]: https://www.r-project.org
[RMarkdown]: http://rmarkdown.rstudio.com
