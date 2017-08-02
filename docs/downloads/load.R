library(bcbioRnaseq)
bcb <- loadRun(
    uploadDir = file.path("data", "final"),
    interestingGroups = c("genotype", "treatment"),
    experimentName = "",
    researcher = "",
    principalInvestigator = "",
    author = getOption("author"),
    email = getOption("email"))
saveData(bcb)
