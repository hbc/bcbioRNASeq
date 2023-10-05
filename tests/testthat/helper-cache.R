lst <- AcidDevTools::cacheTestFiles(
    pkg = .pkgName,
    files = c(
        "bcb_fast.rds",
        "fastrnaseq.tar.gz"
    )
)
cacheDir <- lst[["cacheDir"]]
rm(lst)
