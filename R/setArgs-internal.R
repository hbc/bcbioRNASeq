.setArgs <- function(
    args,
    call,
    fun,
    remove = NULL
) {
    assert_is_list(args)
    # assert_has_names(args)
    # assert_is_call(call)
    # assert_is_function(fun)
    assert_is_any_of(remove, c("character", "NULL"))

    callArgs <- call %>%
        as.list() %>%
        .[-1L] %>%
        .[setdiff(names(.), names(args))]
    args <- c(args, callArgs)

    formalArgs <- fun %>%
        formals() %>%
        .[setdiff(names(.), names(args))]
    args <- c(args, formalArgs)

    args <- args[setdiff(
        x = names(args),
        y = c(remove, "...")
    )]

    stopifnot(!any(duplicated(names(args))))
    args
}
