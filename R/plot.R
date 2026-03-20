#' @export
print.modeller = function(x)
{
    time_var = names(x)[1]
    x2 = tidyr::pivot_longer(x, cols = -tidyr::all_of(time_var))
    x2$name = factor(x2$name, unique(x2$name))

    print(ggplot2::ggplot(x2) +
        ggplot2::geom_line(ggplot2::aes(x = t, y = value, colour = name)))
}
