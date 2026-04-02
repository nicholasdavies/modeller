nice_breaks = scales::breaks_extended(Q = c(1, 5, 2))
# nice_labels = scales::label_number(big.mark = ",",
#     scale_cut = c(" " = 0, " million" = 1e6, " billion" = 1e9, " trillion" = 1e12))

nice_labels = function(x)
{
    dmax = 1e15
    dmin = 1e-3

    resolution = function(x)
    {
        if (length(x)) {
            for (d in 0:5) {
                if (all(abs(round(x, d) - x) < .Machine$double.eps * 100)) {
                    return (10^-d)
                }
            }
        }
        1e-6
    }
    x_valid = x[!is.na(x)]
    res = resolution(x_valid[abs(x_valid) < dmax & abs(x_valid) > dmin & x_valid != 0])

    sc = c(" " = 0, " million" = 1e6, " billion" = 1e9, " trillion" = 1e12)
    names(sc)[1] = ""
    labels = lapply(x, function(x) {
        if (is.na(x)) {
            return ("")
        } else if (x == 0) {
            0
        } else if ((abs(x) < dmax && abs(x) > dmin)) {
            scales::number(x, accuracy = res, big.mark = ",", scale_cut = sc)
        } else {
            exponent = floor(log10(abs(x)))
            mantissa = x / (10^exponent)
            rlang::expr(!!mantissa %*% 10^!!exponent)
        }
    })

    do.call(expression, labels)
}

nice_scales = list(
    ggplot2::scale_x_continuous(breaks = nice_breaks, labels = nice_labels,
        expand = ggplot2::expansion(mult = 1e-6)),
    ggplot2::scale_y_continuous(breaks = nice_breaks, labels = nice_labels,
        expand = ggplot2::expansion(mult = c(1e-3, 0.02)))
)

resolve_palette = function(palette) {
    if (is.function(palette)) return(palette)
    viridis_palettes = c("viridis", "magma", "plasma", "inferno")
    base_palettes = c("Okabe-Ito", "R4", "Tableau 10")
    if (palette %in% viridis_palettes) {
        scales::viridis_pal(option = palette)
    } else if (palette %in% base_palettes) {
        function(n) unname(palette.colors(n, palette = palette))
    } else if (palette == "ggplot2") {
        scales::hue_pal()
    } else {
        scales::brewer_pal(palette = palette)
    }
}

#' Plot model results
#'
#' Plots compartment dynamics from a model simulation. Cumulative counters
#' (compartments with names starting `total_`) are automatically converted
#' to incidence rates. Overlay data points can be added for comparison with
#' observed data.
#'
#' @param x A `model_result` object returned by [run_model()].
#' @param series Character vector of series names to display. `NULL` (the
#'   default) shows all compartments and incidence series.
#' @param overlay Optional data frame of observed data to overlay as points.
#'   Should contain a time column and one or more columns matching
#'   compartment or incidence names.
#' @param palette Colour palette: a palette name (e.g. `"Set1"`, `"Okabe-Ito"`,
#'   `"viridis"`) or a function that takes `n` and returns `n` colours.
#' @param legend Legend position: `"right"`, `"left"`, `"top"`, `"bottom"`,
#'   or `"none"`.
#' @param time_col Name of the time column in `overlay`. Default `"t"`.
#' @param vline Optional x-intercept for a vertical reference line.
#' @param ... Ignored.
#'
#' @return The [ggplot2::ggplot] object (returned by `print()`).
#'
#' @examples
#' m = ode_model(
#'     init = list(S = 999, I = 1, R = 0),
#'     params = list(beta = 0.3, gamma = 0.1),
#'     model = function() list(-beta * S * I, beta * S * I - gamma * I, gamma * I),
#'     times = list(duration = 50)
#' )
#' result = run_model(m)
#' plot(result)
#' plot(result, series = c("I", "R"), palette = "viridis")
#'
#' @export
plot.model_result = function(x, series = NULL, overlay = NULL,
    palette = "Set1", legend = "right", time_col = "t", vline = NULL, ...)
{
    # Separate total_ columns for incidence
    total_cols = grep("^total_", names(x), value = TRUE)
    plot_cols = setdiff(names(x), total_cols)
    main_data = x[plot_cols]

    # Compute incidence from total_ columns
    inc_data = NULL
    if (length(total_cols) > 0) {
        mat = as.matrix(x[c("t", total_cols)])
        grid_step = if (attr(x, "dt") <= 1) 1 else attr(x, "dt")
        grid = seq(min(mat[, 1]), max(mat[, 1]), by = grid_step)
        summ = grid_summarize(mat, grid)
        inc_data = data.frame(t = grid)
        for (col in total_cols) {
            new_name = sub("^total_", "", col)
            vals = summ[, col] / grid_step
            inc_data[[new_name]] = c(vals[1], diff(vals))
        }
    }

    # Determine all series names (used for consistent factor levels / colours)
    compartments = setdiff(names(main_data), "t")
    inc_names = if (!is.null(inc_data)) setdiff(names(inc_data), "t") else character(0)
    all_names = c(compartments, inc_names)

    # Filter to requested series (but keep all_names for colour consistency)
    visible_compartments = if (!is.null(series)) intersect(compartments, series) else compartments
    visible_inc = if (!is.null(series)) intersect(inc_names, series) else inc_names

    # Pivot main data
    if (length(visible_compartments) > 0) {
        plot_data = tidyr::pivot_longer(main_data[c("t", visible_compartments)], cols = -1)
        plot_data$name = factor(plot_data$name, levels = all_names)
    } else {
        plot_data = data.frame(t = numeric(0), name = factor(levels = all_names), value = numeric(0))
    }

    # Choose geom from attribute: "step" for discrete (SSA), "line" for continuous (ODE)
    geom = if (identical(attr(x, "geom"), "step")) ggplot2::geom_step else ggplot2::geom_line

    p = ggplot2::ggplot(plot_data) +
        geom(ggplot2::aes(x = t, y = value, colour = name), linewidth = 1) +
        ggplot2::labs(x = "Time", y = NULL, colour = NULL) +
        nice_scales +
        cowplot::theme_cowplot(font_size = 12) +
        ggplot2::theme(
            axis.line = ggplot2::element_line(linewidth = 0.3),
            axis.ticks = ggplot2::element_line(linewidth = 0.3),
            legend.position = legend
        )

    # Vertical line
    if (!is.null(vline)) {
        p = p + ggplot2::geom_vline(xintercept = vline)
    }

    # Overlay incidence
    if (!is.null(inc_data) && length(visible_inc) > 0) {
        inc_long = tidyr::pivot_longer(inc_data[c("t", visible_inc)], cols = -1)
        inc_long$name = factor(inc_long$name, levels = all_names)
        if (nrow(inc_long) > 0) {
            p = p + geom(data = inc_long,
                ggplot2::aes(x = t, y = value, colour = name), linewidth = 1)
        }
    }

    # Overlay data points
    if (!is.null(overlay)) {
        other_cols = setdiff(names(overlay), time_col)
        overlay_long = tidyr::pivot_longer(overlay, cols = dplyr::all_of(other_cols),
            names_to = "name", values_to = "value")
        overlay_long$name = factor(overlay_long$name, levels = all_names)
        if (!is.null(series)) {
            overlay_long = overlay_long[overlay_long$name %in% series, ]
        }
        p = p + ggplot2::geom_point(data = overlay_long,
            ggplot2::aes(x = .data[[time_col]], y = value, colour = name),
            size = 2.5) +
            ggplot2::geom_point(data = overlay_long,
            ggplot2::aes(x = .data[[time_col]], y = value),
            shape = 1, colour = "black", size = 2.5, stroke = 0.4)
    }

    # Apply palette
    palette_fn = resolve_palette(palette)
    n_colours = length(all_names)
    p = p + ggplot2::scale_colour_manual(values = palette_fn(n_colours))

    p
}

#' @export
print.model_result = function(x, ...) {
    compartments = setdiff(names(x), "t")
    t_range = range(x$t)
    cat(sprintf("Model result: %d rows, t = [%g, %g], compartments: %s\n",
        nrow(x), t_range[1], t_range[2], paste(compartments, collapse = ", ")))
    print.data.frame(x, ...)
}
