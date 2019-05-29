# Contributing to development

Please follow all [recommended guidelines](https://basejump.acidgenomics.com/basejump/CONTRIBUTING.html) defined in the [basejump](https://basejump.acidgenomics.com) package.

## Administrator guidelines

- Do not commit changes directly to master branch.
- Please avoid creating unnecessary new branches. Use a fork instead.

## R Markdown templates

Check that all templates render correctly.

```r
git clone https://github.com/bcbio/bcbio_rnaseq_output_example.git
cd bcbio_rnaseq_output_example
Rscript -e 'testthat::test_file("test_reports.R")'
```

Ensure that all expected HTML files are created inside the `reports/` directory, and that all figures and tables look correct.
