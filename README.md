# bcbioRNASeq

[![Travis CI](https://travis-ci.org/hbc/bcbioRNASeq.svg?branch=master)](https://travis-ci.org/hbc/bcbioRNASeq)
[![AppVeyor CI](https://ci.appveyor.com/api/projects/status/s0rrc28fwr0ua2wr/branch/master?svg=true)](https://ci.appveyor.com/project/mjsteinbaugh/bcbiornaseq/branch/master)
[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/r-bcbiornaseq/badges/version.svg)](https://anaconda.org/bioconda/r-bcbiornaseq)

[R][] package for [bcbio][] RNA-seq analysis.

## Workflow paper

Steinbaugh MJ, Pantano L, Kirchner RD, Barrera V, Chapman BA, Piper ME, Mistry M, Khetani RS, Rutherford KD, Hoffman O, Hutchinson JN, Ho Sui SJ. (2018). [bcbioRNASeq: R package for bcbio RNA-seq analysis.][workflow paper] *F1000Research* 6:1976.

```r
citation("bcbioRNASeq")
```

## Installation

### [Bioconductor][] method

We recommend installing the package with [BiocManager][].

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
library(BiocManager)
install("remotes")
install("hbc/bcbioRNASeq")
```

### [conda][] method

Configure [conda][] to use the [bioconda][] channels.

```sh
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

```sh
conda create --name r
conda activate r
```

Launch [R][] and check that it is set up correctly with the `capabilities()` function. Note that `X11 = TRUE` is required for graphical output, and requires X11 forwarding over SSH.

Now you're ready to install `r-bcbiornaseq`.

```sh
conda install -c bioconda r-bcbiornaseq
```

Note that there is currently a bug with [conda][] and `libgfortran`. You may need to install `libgfortran-ng` to get the bcbioRNASeq package to load in [R][].

```sh
conda install libgfortran-ng
```

## Load [bcbio][] RNA-seq data

```r
library(bcbioRNASeq)
bcb <- bcbioRNASeq(
    uploadDir = "bcbio_rnaseq_run/final",
    interestingGroups = c("genotype", "treatment"),
    organism = "Homo sapiens"
)
saveData(bcb, dir = ".")
```

This will return a `bcbioRNASeq` object, which is an extension of the [Bioconductor][] [RangedSummarizedExperiment][] container class. Consult the `bcbioRNASeq()` constructor function documentation for detailed information on the supported parameters:

```r
help(topic = "bcbioRNASeq", package = "bcbioRNASeq")
```

### Sample metadata

When loading a [bcbio][] RNA-seq run, the sample metadata will be imported automatically from the `project-summary.yaml` file in the final upload directory. If you notice any typos in your metadata after completing the run, these can be corrected by editing the YAML file. Alternatively, you can pass in a sample metadata file into `bcbioRNASeq()` using the `sampleMetadataFile` argument.

#### Metadata file example

The samples in the [bcbio][] run must map to the `description` column. The values provided in `description` must be unique. These values will be sanitized into syntactically valid names (see `help("make.names")`), and assigned as the column names of the `bcbioRNASeq` object. The original values are stored as the `sampleName` column in `colData()`, and are used for all plotting functions.

| description | genotype |
|-------------|----------|
| sample1     | wildtype |
| sample2     | knockout |
| sample3     | wildtype |
| sample4     | knockout |

## [R Markdown][] templates

The package provides multiple [R Markdown][] templates, including quality control, differential expression using [DESeq2][], and functional enrichment analysis. These are available in [RStudio][] at `File` -> `New File` -> `R Markdown...` -> `From Template`.

## References

The papers and software cited in our workflows are available as a [shared library](https://paperpile.com/shared/e1q8fn) on [Paperpile][].

[bcbio]: https://github.com/chapmanb/bcbio-nextgen
[BiocManager]: https://cran.r-project.org/package=BiocManager
[bioconda]: https://bioconda.github.io
[Bioconductor]: https://bioconductor.org
[conda]: https://conda.io
[DESeq2]: https://doi.org/doi:10.18129/B9.bioc.DESeq2
[Paperpile]: https://paperpile.com
[R]: https://www.r-project.org
[R Markdown]: http://rmarkdown.rstudio.com
[RStudio]: https://www.rstudio.com
[RangedSummarizedExperiment]: https://doi.org/doi:10.18129/B9.bioc.SummarizedExperiment
[Workflow paper]: https://f1000research.com/articles/6-1976/v2
