## Testing both gene and transcript level. GFF3 files are also supported, but
## we're only testing GTF here for speed. This functionality is covered in
## basejump tests also.
test_that("GFF/GTF file", {
    for (level in eval(formals(bcbioRNASeq)[["level"]])) {
        gffURL <- pasteUrl(
            "ftp.ensembl.org",
            "pub",
            "release-90",
            "gtf",
            "mus_musculus",
            "Mus_musculus.GRCm38.90.gtf.gz",
            protocol = "ftp"
        )
        gffFile <- file.path(cacheDir, basename(gffURL))
        if (!file.exists(gffFile)) {
            initDir(cacheDir)
            download.file(url = gffURL, destfile = gffFile)
        }
        object <- bcbioRNASeq(
            uploadDir = uploadDir,
            level = level,
            organism = organism,
            gffFile = gffFile
        )
        expect_s4_class(object, "bcbioRNASeq")
    }
})
