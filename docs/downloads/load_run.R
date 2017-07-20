library(bcbioRnaseq)
bcb <- load_run(
    upload_dir = file.path("data", "final"),
    sample_metadata_file = file.path("meta", "sample_metadata.xlsx"),
    interesting_groups = c("genotype", "treatment"),
    experiment_name = "",
    researcher = "",
    principal_investigator = "",
    author = getOption("author"),
    email = getOption("email"))
save_data(bcb)
