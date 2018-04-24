# bcbioRNASeq

[![Travis CI](https://travis-ci.org/hbc/bcbioRNASeq.svg?branch=master)](https://travis-ci.org/hbc/bcbioRNASeq)
[![AppVeyor CI](https://ci.appveyor.com/api/projects/status/s0rrc28fwr0ua2wr/branch/master?svg=true)](https://ci.appveyor.com/project/mjsteinbaugh/bcbiornaseq/branch/master)
[![Codecov](https://codecov.io/gh/hbc/bcbioRNASeq/branch/master/graph/badge.svg)](https://codecov.io/gh/hbc/bcbioRNASeq)
[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/r-bcbiornaseq/badges/version.svg)](https://anaconda.org/bioconda/r-bcbiornaseq)

Quality control and differential expression for [bcbio][] RNA-seq experiments.


## Installation

This is an [R][] package.

### [Bioconductor][] method

```r
source("https://bioconductor.org/biocLite.R")
biocLite("devtools")
biocLite(
    "hbc/bcbioRNASeq",
    dependencies = c("Depends", "Imports", "Suggests")
)
```

### [conda][] method

```bash
conda install -c bioconda r-bcbiornaseq
```


## Load [bcbio][] run

```r
library(bcbioRNASeq)
bcb <- bcbioRNASeq(
    uploadDir = "bcbio_rnaseq_run/final",
    interestingGroups = c("genotype", "treatment"),
    organism = "Homo sapiens"
)
# Back up all data inside bcbioRNASeq object
flatFiles <- flatFiles(bcb)
saveData(bcb, flatFiles)
```

This will return a `bcbioRNASeq` object, which is an extension of the [Bioconductor][] [RangedSummarizedExperiment][RSE] container class.

Parameters:

- `uploadDir`: Path to the [bcbio][] final upload directory.
- `interestingGroups`: Character vector of the column names of interest in the sample metadata, which is stored in the `colData()` accessor slot of the `bcbioRNASeq` object. These values should be formatted in camelCase, and can be reassigned in the object after creation (e.g. `interestingGroups(bcb) <- c("batch", "age")`). They are used for data visualization in the quality control utility functions.
- `organism`: Organism name. Use the full latin name (e.g. "Homo sapiens").

Consult `help("bcbioRNASeq", "bcbioRNASeq")` for additional documentation.


## [R Markdown][] templates

This package provides multiple [R Markdown][] templates, including Quality Control and Differential Expression using [DESeq2][], which are available in [RStudio][] at `File` -> `New File` -> `R Markdown...` -> `From Template`.

### Examples

View example [HTML reports](http://bcb.io/bcbio_rnaseq_output_example) rendered from the default [R Markdown][] templates included in the package:

- [Quality Control](http://bcb.io/bcbio_rnaseq_output_example/qc-master.html)
- [Differential Expression](http://bcb.io/bcbio_rnaseq_output_example/de-master.html)


## Sample metadata

When loading a [bcbio][] RNA-seq run, the sample metadata will be imported automatically from the `project-summary.yaml` file in the final upload directory. If you notice any typos in your metadata after completing the run, these can be corrected by editing the YAML file. Alternatively, you can pass in a sample metadata file into `bcbioRNASeq()` using the `sampleMetadataFile` argument.

### Metadata file example

The samples in the [bcbio][] run must map to the `description` column. The values provided in the `description` column must be unique, and can contain any characters, including spaces or hyphens. These values will be sanitized into syntactically valid names (see `help("makeNames", "basejump")` for more information), and assigned as the column names (`colnames()`) of the `bcbioRNASeq` object. The unmodified values are accessible from the `sampleName` column in `colData()`.

| description | genotype |
|-------------|----------|
| sample1     | wildtype |
| sample2     | knockout |
| sample3     | wildtype |
| sample4     | knockout |


## Citation

```r
citation("bcbioRNASeq")
```

Steinbaugh MJ, Pantano L, Kirchner RD, Barrera V, Chapman BA, Piper ME, Mistry M, Khetani RS, Rutherford KD, Hoffman O, Hutchinson JN, Ho Sui SJ. (2017). [bcbioRNASeq: R package for bcbio RNA-seq analysis.](https://f1000research.com/articles/6-1976/v1) *F1000Research* 6:1976.


## References

The papers and software cited in our workflows are available as a [shared library](https://paperpile.com/shared/e1q8fn) on [Paperpile][].


[bcbio]: https://github.com/chapmanb/bcbio-nextgen
[Bioconductor]: https://bioconductor.org
[conda]: https://conda.io
[DESeq2]: https://doi.org/doi:10.18129/B9.bioc.DESeq2
[Paperpile]: https://paperpile.com
[R]: https://www.r-project.org
[R Markdown]: http://rmarkdown.rstudio.com
[RStudio]: https://www.rstudio.com
[RSE]: https://doi.org/doi:10.18129/B9.bioc.SummarizedExperiment
