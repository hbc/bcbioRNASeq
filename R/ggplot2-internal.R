# nolint start
#
# Reverse the order of a categorical axis in ggplot2.
# - https://gist.github.com/jennybc/6f3fa527b915b920fdd5
# - https://gist.github.com/jennybc/6f3fa527b915b920fdd5#gistcomment-2709675
#
# Useful alternate approach:
# > aes(x = reorder(the_factor, desc(the_factor)), ...)
#
# nolint end

.flipMode <- function(object) {
    assert(is(object, "ggplot"))
    data <- object[["data"]]
    assert(is.data.frame(data))
    samples <- data[["sampleName"]]
    assert(is.factor(samples))
    object +
        scale_x_discrete(limits = rev(levels(samples))) +
        coord_flip()
}
