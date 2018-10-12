data <- tibble(
    description = paste0("sample", seq_len(4L)),
    genotype = rep(c("wildtype", "knockout"), times = 2L)
)
bb8::tabular(data)
