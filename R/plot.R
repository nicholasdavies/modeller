nice_breaks = scales::breaks_extended(Q = c(1, 5, 2))
nice_labels = scales::label_number(big.mark = ",",
    scale_cut = c(" " = 0, " million" = 1e6, " billion" = 1e9, " trillion" = 1e12))

nice_scales = list(
    ggplot2::scale_x_continuous(breaks = nice_breaks, labels = nice_labels,
        expand = ggplot2::expansion(add = 0)),
    ggplot2::scale_y_continuous(breaks = nice_breaks, labels = nice_labels,
        expand = ggplot2::expansion(mult = c(0, 0.02)))
)
