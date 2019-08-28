# bcbioRNASeq

[![Repo status: active](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Travis CI build status](https://travis-ci.org/hbc/bcbioRNASeq.svg?branch=master)](https://travis-ci.org/hbc/bcbioRNASeq)
[![AppVeyor CI build status](https://ci.appveyor.com/api/projects/status/s0rrc28fwr0ua2wr/branch/master?svg=true)](https://ci.appveyor.com/project/mjsteinbaugh/bcbiornaseq/branch/master)
[![Anaconda version](https://anaconda.org/bioconda/r-bcbiornaseq/badges/version.svg) ![Anaconda latest release date](https://anaconda.org/bioconda/r-bcbiornaseq/badges/latest_release_date.svg) ![Anaconda downloads](https://anaconda.org/bioconda/r-bcbiornaseq/badges/downloads.svg)](https://anaconda.org/bioconda/r-bcbiornaseq)

[R][] package for [bcbio][] RNA-seq analysis.

## Workflow paper

Steinbaugh MJ, Pantano L, Kirchner RD, Barrera V, Chapman BA, Piper ME, Mistry M, Khetani RS, Rutherford KD, Hoffman O, Hutchinson JN, Ho Sui SJ. (2018). [bcbioRNASeq: R package for bcbio RNA-seq analysis.][workflow paper] *F1000Research* 6:1976.

```r
citation("bcbioRNASeq")
```

## Installation

### [R][] method

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
}
Sys.setenv(R_REMOTES_UPGRADE = "always")
# Set `GITHUB_PAT` in `~/.Renviron` if you get a rate limit error.
remotes::install_github("hbc/bcbioRNASeq")
```

Here's how to update to the latest version on GitHub:

```r
Sys.setenv(R_REMOTES_UPGRADE = "always")
remotes::update_packages()
```

Always check that your Bioconductor installation is valid before proceeding.

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::valid()
```

### [conda][] method

Configure [conda][] to use the [bioconda][] channels.

```sh
conda install -c bioconda r-bcbiornaseq
```

## Load [bcbio][] RNA-seq data

```r
library(bcbioRNASeq)
object <- bcbioRNASeq(
    uploadDir = file.path("bcbio", "final"),
    interestingGroups = c("genotype", "treatment"),
    organism = "Homo sapiens"
)
saveData(object, dir = ".")
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

## Troubleshooting

### Invalid object

If you encounter a `validObject` error when attempting to load a `bcbioRNASeq` object from a previous analysis, run this step to update the object to the current version of the package:

```r
object <- updateObject(object)
validObject(object)
## [1] TRUE
```

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
