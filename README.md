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
biocLite("ensembldb")
biocLite(
    "hbc/bcbioRNASeq",
    dependencies = c("Depends", "Imports", "Suggests")
)
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
- `interestingGroups`: Character vector of the column names of interest in the sample metadata, which is stored in the `colData()` accessor slot of the `bcbioRNASeq` object. These values should be formatted in camelCase, and can be reassigned in the object after creation (e.g. `interestingGroups(bcb) <- c("batch", "age")`). They are used for data visualization in the quality control utility functions.

Consult `help("loadRNASeq", "bcbioRNASeq")` for additional documentation.


## [R Markdown][] templates

This package provides multiple [R Markdown][] templates, including Quality Control and Differential Expression using [DESeq2][], which are available in [RStudio][] at `File` -> `New File` -> `R Markdown...` -> `From Template`.

### Examples

View example [HTML reports](http://bcb.io/bcbio_rnaseq_output_example) rendered from the default [R Markdown][] templates included in the package:

- [Quality Control](http://bcb.io/bcbio_rnaseq_output_example/qc-master.html)
- [Differential Expression](http://bcb.io/bcbio_rnaseq_output_example/de-master.html)


## Citation

Steinbaugh MJ, Pantano L, Kirchner RD, Barrera V, Chapman BA, Piper ME, Mistry M, Khetani RS, Rutherford KD, Hoffman O, Hutchinson JN, Ho Sui SJ. (2017). [bcbioRNASeq: R package for bcbio RNA-seq analysis.](https://f1000research.com/articles/6-1976/v1) *F1000 Research* 6:1976.


## References

The papers and software cited in our workflows are available as a [shared library](https://paperpile.com/shared/e1q8fn) on [Paperpile][].


[bcbio]: https://github.com/chapmanb/bcbio-nextgen
[Bioconductor]: https://bioconductor.org
[DESeq2]: https://bioconductor.org/packages/release/bioc/html/DESeq2.html
[devtools]: https://cran.r-project.org/package=devtools
[Paperpile]: https://paperpile.com
[R]: https://www.r-project.org
[R Markdown]: http://rmarkdown.rstudio.com
[RStudio]: https://www.rstudio.com
[SummarizedExperiment]: http://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html
