Sys.setenv(TZ = "America/New_York")
sessioninfo::session_info()
rcmdcheck::rcmdcheck(
    path = ".",
    args = c(
        "--no-build-vignettes",
        "--no-manual",
        "--no-vignettes",
        "--timings"
    ),
    build_args = c(
        "--no-build-vignettes",
        "--no-manual"
    ),
    error_on = "error"
)

if (packageVersion("base") >= "3.6") {
    BiocCheck::BiocCheck(
        package = ".",
        `quit-with-status` = TRUE
    )
    lintr::lint_package()
}
