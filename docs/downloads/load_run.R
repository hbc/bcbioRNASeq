library(bcbioRnaseq)
bcb <- load_run(
    upload_dir = file.path("data", "final"),
    interesting_groups = c("genotype", "treatment"),
    experiment_name = "",
    researcher = "",
    principal_investigator = "",
    author = getOption("author"),
    email = getOption("email"))
save_data(bcb)
