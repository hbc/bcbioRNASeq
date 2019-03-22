# bcbioRNASeq

[![Travis CI](https://travis-ci.org/hbc/bcbioRNASeq.svg?branch=master)](https://travis-ci.org/hbc/bcbioRNASeq)
[![AppVeyor CI](https://ci.appveyor.com/api/projects/status/s0rrc28fwr0ua2wr/branch/master?svg=true)](https://ci.appveyor.com/project/mjsteinbaugh/bcbiornaseq/branch/master)
[![Codecov](https://codecov.io/gh/hbc/bcbioRNASeq/branch/master/graph/badge.svg)](https://codecov.io/gh/hbc/bcbioRNASeq)
[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/r-bcbiornaseq/badges/version.svg)](https://anaconda.org/bioconda/r-bcbiornaseq)

[R][] package for [bcbio][] RNA-seq analysis.

## Workflow paper

Steinbaugh MJ, Pantano L, Kirchner RD, Barrera V, Chapman BA, Piper ME, Mistry M, Khetani RS, Rutherford KD, Hoffman O, Hutchinson JN, Ho Sui SJ. (2018). [bcbioRNASeq: R package for bcbio RNA-seq analysis.][F1000 paper] *F1000Research* 6:1976.

```r
citation("bcbioRNASeq")
```

## Installation

### [Bioconductor][] method

We recommend installing the package with [BiocManager][].

```r
if (!require("BiocManager")) {
    install.packages("BiocManager")
}
BiocManager::install("remotes")
BiocManager::install("hbc/bcbioRNASeq", ref = "bioconductor")
```

For [R][] < 3.5, [BiocManager][] is not supported. Use `BiocInstaller::biocLite()` instead of `BiocManager::install()`. This requires sourcing the legacy [Bioconductor][] `biocLite.R` script.

```r
# Try `http://` if `https://` URLs are not supported.
source("https://bioconductor.org/biocLite.R")
```

### [conda][] method

Configure [conda][] to use the [bioconda][] channels.

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

To avoid version issues, your `.condarc` file should only contain the following channels, in this order:

```
channels:
  - conda-forge
  - bioconda
  - defaults
```

We recommend installing into a clean [conda][] environment:

```bash
conda create --name r
conda activate r
```

Launch [R][] and check that it is set up correctly with the `capabilities()` function. Note that `X11 = TRUE` is required for graphical output, and requires X11 forwarding over SSH.

Now you're ready to install `r-bcbiornaseq`.

```bash
conda install -c bioconda r-bcbiornaseq
```

Note that there is currently a bug with [conda][] and `libgfortran`. You may need to install `libgfortran-ng` to get the bcbioRNASeq package to load in [R][].

```bash
conda install libgfortran-ng
```

### Pinned [workflow paper][F1000] version

Installing the specific version of [bcbioRNASeq][] used in the [F1000 workflow paper][F1000] is easiest using [conda][].

```bash
conda install r-bcbiornaseq=0.2.4
```

Alternatively, you can use [devtools][] to install the pinned version [bcbioRNASeq][] and its dependencies. First install [bcbioRNASeq][] using [BiocManager][] (see above), and then run these commands.

```r
devtools::install("hbc/bcbioRNASeq", ref = "v0.2.4", dependencies = FALSE)
devtools::install("hbc/bcbioBase", ref = "v0.2.12", dependencies = FALSE)
devtools::install("steinbaugh/basejump", ref = "v0.5.0", dependencies = FALSE)
```

## Loading [bcbio][] RNA-seq data

```r
library(bcbioRNASeq)
bcb <- bcbioRNASeq(uploadDir = "bcbio/final")
saveData(bcb, dir = ".")
```

This will return a `bcbioRNASeq` object, which is an extension of the [Bioconductor][] [RangedSummarizedExperiment][] container class. Consult the `bcbioRNASeq()` generator function documentation for detailed information on the supported parameters.

```r
help(topic = "bcbioRNASeq", package = "bcbioRNASeq")
```

## [R Markdown][] templates

The package provides multiple [R Markdown][] templates, including quality control, differential expression using [DESeq2][], and functional enrichment analysis. They are available in [RStudio][] at `File` -> `New File` -> `R Markdown...` -> `From Template`.

### Example reports

- [Quality Control](http://bcb.io/bcbio_rnaseq_output_example/qc-master.html)
- [Differential Expression](http://bcb.io/bcbio_rnaseq_output_example/de-master.html)
- [Functional Analysis](http://bcb.io/bcbio_rnaseq_output_example/fa-master.html)

## References

The papers and software cited in our workflows are available as a [shared library](https://paperpile.com/shared/e1q8fn) on [Paperpile][].

[bcbio]: https://bcbio-nextgen.readthedocs.io/
[bcbioRNASeq]: https://bioinformatics.sph.harvard.edu/bcbioRNASeq/
[BiocManager]: https://cran.r-project.org/package=BiocManager
[bioconda]: https://bioconda.github.io/
[Bioconductor]: https://bioconductor.org/
[conda]: https://conda.io/
[DESeq2]: https://bioconductor.org/packages/DESeq2/
[devtools]: https://cran.r-project.org/package=devtools
[F1000 paper]: http://dx.doi.org/10.12688/f1000research.12093.2
[Paperpile]: https://paperpile.com/
[R]: https://www.r-project.org/
[R Markdown]: https://rmarkdown.rstudio.com/
[RStudio]: https://www.rstudio.com/
[RangedSummarizedExperiment]: https://bioconductor.org/packages/SummarizedExperiment/
