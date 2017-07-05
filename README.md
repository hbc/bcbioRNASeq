[Bioconductor]: https://bioconductor.org
[bcbio-nextgen]: https://github.com/chapmanb/bcbio-nextgen
[devtools]: https://cran.r-project.org/package=devtools
[R]: https://www.r-project.org



# bcbioRnaseq

[![Build Status](https://travis-ci.org/hbc/bcbioRnaseq.svg?branch=master)](https://travis-ci.org/hbc/bcbioRnaseq)
[![Project Status: WIP - Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)
[![codecov](https://codecov.io/gh/hbc/bcbioRnaseq/branch/master/graph/badge.svg)](https://codecov.io/gh/hbc/bcbioRnaseq)

Quality control and differential expression for [bcbio-nextgen][] RNA-seq experiments.


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

[HTML reports](http://bcb.io/bcbio_rnaseq_output_example) rendered from the
default RMarkdown templates included in the package.

* [QC example](http://bcb.io/bcbio_rnaseq_output_example/qc.html)
* [DE example](http://bcb.io/bcbio_rnaseq_output_example/de.html)


## Contribute

- Only current release branch will be maintened 
- For major changes, create a new branch and a pull request that will be revised.
- For small changes, commit directly to master.
- For fixes in releases, commit to releases branches.
- Open issues or make pull request for changes affecting `bcbioRNADataSet` object or `setup.R`
