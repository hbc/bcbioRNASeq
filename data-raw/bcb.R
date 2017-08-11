library(bcbioRnaseq)
path <- system.file("extra", package = "bcbioRnaseq")
bcb <- loadRun(file.path(path, "bcbio"))
saveData(bcb, compress = "xz")
