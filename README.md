# bcbioRNASeq

[![Install with Bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](https://bioconda.github.io/recipes/r-bcbiornaseq/README.html)
![Lifecycle: retired](https://img.shields.io/badge/lifecycle-retired-red.svg)

**NOTE: [bcbio-nextgen][bcbio] is no longer under active development.**
Refer to the [notice of discontinuation][] for additional details.

[R][] package for [bcbio][] RNA-seq analysis.

## Workflow paper

Steinbaugh MJ, Pantano L, Kirchner RD, Barrera V, Chapman BA, Piper ME, Mistry
M, Khetani RS, Rutherford KD, Hoffman O, Hutchinson JN, Ho Sui SJ. (2018).
[bcbioRNASeq: R package for bcbio RNA-seq analysis.][workflow paper]
_F1000Research_ 6:1976.

```r
citation("bcbioRNASeq")
```

## Installation

This is an [R][] package.

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
install.packages(
    pkgs = "bcbioRNASeq",
    repos = c(
        "https://r.acidgenomics.com",
        BiocManager::repositories()
    ),
    dependencies = TRUE
)
```

### [Conda][] method

Configure [Conda][] to use the [Bioconda][] channels.

```sh
# Don't install recipe into base environment.
conda create --name='r-bcbiornaseq' 'r-bcbiornaseq'
conda activate 'r-bcbiornaseq'
R
```

## Load [bcbio][] RNA-seq data

```r
library(bcbioRNASeq)
object <- bcbioRNASeq(
    uploadDir = file.path("bcbio", "final"),
    interestingGroups = c("genotype", "treatment"),
    organism = "Homo sapiens"
)
```

This will return a `bcbioRNASeq` object, which is an extension of the
[Bioconductor][] [RangedSummarizedExperiment][] container class. Consult the
`bcbioRNASeq()` constructor function documentation for detailed information on
the supported parameters:

```r
help(topic = "bcbioRNASeq", package = "bcbioRNASeq")
```

### Sample metadata

When loading a [bcbio][] RNA-seq run, the sample metadata will be imported
automatically from the `project-summary.yaml` file in the final upload
directory. If you notice any typos in your metadata after completing the run,
these can be corrected by editing the YAML file. Alternatively, you can pass in
a sample metadata file into `bcbioRNASeq()` using the `sampleMetadataFile`
argument.

#### Metadata file example

The samples in the [bcbio][] run must map to the `description` column. The
values provided in `description` must be unique. These values will be sanitized
into syntactically valid names (see `help("make.names")`), and assigned as the
column names of the `bcbioRNASeq` object. The original values are stored as the
`sampleName` column in `colData()`, and are used for all plotting functions.

| description | genotype |
| ----------- | -------- |
| sample1     | wildtype |
| sample2     | knockout |
| sample3     | wildtype |
| sample4     | knockout |

## Differential expression

We've designed bcbioRNASeq to easily hand off to [DESeq2][], [edgeR][], or
[limma-voom][] for differential expression analysis.

DESeq2: Coerce `bcbioRNASeq` to `DESeqDataSet`.

```r
dds <- as.DESeqDataSet(object)
```

edgeR or limma-voom: Coerce `bcbioRNASeq` to `DGEList`.

```r
dge <- as.DGEList(object)
```

## [R Markdown][] templates

The package provides multiple [R Markdown][] templates, including quality
control, differential expression using [DESeq2][], and functional enrichment
analysis. These are available in [RStudio][] at `File` -> `New File` ->
`R Markdown...` -> `From Template`.

## Troubleshooting

### Invalid object

If you encounter a `validObject` error when attempting to load a `bcbioRNASeq`
object from a previous analysis, run this step to update the object to the
current version of the package:

```r
object <- updateObject(object)
validObject(object)
## [1] TRUE
```

## References

The papers and software cited in our workflows are available as a [shared
library](https://paperpile.com/shared/e1q8fn) on [Paperpile][].

[bcbio]: https://github.com/bcbio/bcbio-nextgen/
[bioconda]: https://bioconda.github.io/
[bioconductor]: https://bioconductor.org/
[conda]: https://conda.io/
[deseq2]: https://bioconductor.org/packages/DESeq2/
[edger]: https://bioconductor.org/packages/edgeR/
[limma-voom]: https://bioconductor.org/packages/limma/
[notice of discontinuation]: https://github.com/bcbio/bcbio-nextgen/issues/3749
[paperpile]: https://paperpile.com/
[r markdown]: https://rmarkdown.rstudio.com/
[r]: https://www.r-project.org/
[rangedsummarizedexperiment]: https://bioconductor.org/packages/SummarizedExperiment/
[rstudio]: https://posit.co/products/open-source/rstudio/
[workflow paper]: https://doi.org/10.12688/f1000research.12093.2
