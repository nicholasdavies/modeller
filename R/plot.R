# Wrap the closure so labeling::extended is only referenced at runtime,
# not embedded in the package namespace (which would require declaring
# labeling in Imports).
nice_breaks = function(x) scales::breaks_extended(Q = c(1, 5, 2))(x)
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
        expand = ggplot2::expansion(mult = 1e-3)),
    ggplot2::scale_y_continuous(breaks = nice_breaks, labels = nice_labels,
        expand = ggplot2::expansion(mult = c(1e-3, 0.02)))
)

# Extend a fixed set of colours to exactly n. Within range the first n are
# returned unchanged; beyond the palette's natural maximum, interpolate so
# series stay distinct rather than recycling or padding with NA.
extend_colours = function(cols, n) {
    cols = cols[!is.na(cols)]
    if (n <= length(cols)) cols[seq_len(n)] else grDevices::colorRampPalette(cols)(n)
}

resolve_palette = function(palette) {
    if (is.function(palette)) return(palette)
    viridis_palettes = c("viridis", "magma", "plasma", "inferno")
    base_palettes = c("Okabe-Ito", "R4", "Tableau 10")
    if (palette == "Grey") {
        function(n) rep(c("#000000", "#7F7F7F"), length.out = n)
    } else if (palette %in% viridis_palettes) {
        scales::viridis_pal(option = palette)
    } else if (palette %in% base_palettes) {
        # palette.colors() with no n returns the palette's full colour set.
        full = unname(grDevices::palette.colors(palette = palette))
        function(n) extend_colours(full, n)
    } else if (palette == "ggplot2") {
        scales::hue_pal()
    } else {
        # ColorBrewer palettes cap at brewer.pal.info$maxcolors; build the full
        # set once, then interpolate past it.
        max_n = RColorBrewer::brewer.pal.info[palette, "maxcolors"]
        full = scales::brewer_pal(palette = palette)(max_n)
        function(n) extend_colours(full, n)
    }
}

# Linetypes for the "Grey" palette: paired with the alternating colours so
# that each pair of black + grey shares a linetype, then the linetype
# advances. Solid, dashed, dotted, dotdash; recycled if more series.
grey_linetypes = function(n) {
    rep(c("solid", "solid", "dashed", "dashed",
          "dotted", "dotted", "dotdash", "dotdash"), length.out = n)
}

#' Plot model results
#'
#' Plots compartment dynamics from a model simulation. Overlay data points
#' can be added for comparison with observed data.
#'
#' @param x A `model_result` object returned by [run_model()].
#' @param series Character vector of series names to display. `NULL` (the
#'   default) shows all compartments and incidence series.
#' @param overlay Optional data frame of observed data to overlay as points.
#'   Should contain a time column and one or more columns matching
#'   compartment or incidence names.
#' @param palette Colour palette: a palette name (e.g. `"Dark2"`, `"Okabe-Ito"`,
#'   `"viridis"`) or a function that takes `n` and returns `n` colours.
#' @param legend Legend position: `"right"`, `"left"`, `"top"`, `"bottom"`,
#'   or `"none"`.
#' @param time_col Name of the time column in `overlay`. Default `"t"`.
#' @param vline Optional x-intercept for a vertical reference line.
#' @param title Plot title. `NULL` (default) for no title.
#' @param xlab X-axis label. Default `"Time"`. `NULL` for no label.
#' @param ylab Y-axis label. `NULL` (default) for no label.
#' @param compare Optional second `model_result` to draw underneath in
#'   thinner, fainter lines for visual comparison.
#' @param ... Ignored.
#'
#' @return The [ggplot2::ggplot] object (returned by `print()`).
#'
#' @examples
#' m = ode_model(
#'     init = list(S = 999, I = 1, R = 0),
#'     params = list(beta = 0.3, gamma = 0.1, time = 50),
#'     equations = function() {
#'         N = S + I + R
#'         d(S) = -beta * I / N * S
#'         d(I) =  beta * I / N * S - gamma * I
#'         d(R) =  gamma * I
#'     }
#' )
#' result = run_model(m)
#' plot(result)
#' plot(result, series = c("I", "R"), palette = "viridis")
#'
#' @export
plot.model_result = function(x, series = NULL, overlay = NULL,
    palette = "Dark2", legend = "right", time_col = "t", vline = NULL,
    title = NULL, xlab = "Time", ylab = NULL, compare = NULL, ...)
{
    # Determine all series names (used for consistent factor levels / colours)
    all_names = setdiff(names(x), "t")
    visible = if (!is.null(series)) intersect(all_names, series) else all_names

    # Pivot to long format
    if (length(visible) > 0) {
        plot_data = tidyr::pivot_longer(x[c("t", visible)], cols = -1)
        plot_data$name = factor(plot_data$name, levels = all_names)
    } else {
        plot_data = data.frame(t = numeric(0), name = factor(levels = all_names), value = numeric(0))
    }

    # Choose geom from attribute: "step" for discrete (SSA), "line" for continuous (ODE)
    geom = if (identical(attr(x, "geom"), "step")) ggplot2::geom_step else ggplot2::geom_line

    # The "Grey" palette also varies linetype with series so series remain
    # distinguishable in monochrome.
    use_linetypes = identical(palette, "Grey")
    line_aes = if (use_linetypes) {
        ggplot2::aes(x = t, y = value, colour = name, linetype = name)
    } else {
        ggplot2::aes(x = t, y = value, colour = name)
    }

    p = ggplot2::ggplot() +
        ggplot2::geom_blank(data = data.frame(t = x$t[1], value = 0),
            ggplot2::aes(x = t, y = value)) +
        ggplot2::labs(title = title, x = xlab, y = ylab, colour = NULL, linetype = NULL) +
        nice_scales +
        cowplot::theme_cowplot(font_size = 12) +
        ggplot2::theme(
            axis.line = ggplot2::element_line(linewidth = 0.3),
            axis.ticks = ggplot2::element_line(linewidth = 0.3),
            legend.position = legend,
            plot.margin = ggplot2::margin(t = 5, r = 15, b = 5, l = 5)
        )

    # Comparison series, drawn underneath in thinner, fainter lines
    if (!is.null(compare)) {
        compare_visible = intersect(visible, setdiff(names(compare), "t"))
        if (length(compare_visible) > 0) {
            compare_long = tidyr::pivot_longer(compare[c("t", compare_visible)], cols = -1)
            compare_long$name = factor(compare_long$name, levels = all_names)
            compare_geom = if (identical(attr(compare, "geom"), "step")) ggplot2::geom_step else ggplot2::geom_line
            p = p + compare_geom(data = compare_long, line_aes,
                linewidth = 0.6, alpha = 0.7)
        }
    }

    p = p + geom(data = plot_data, line_aes, linewidth = 1)

    # Vertical line with crosses at compartment values
    if (!is.null(vline)) {
        p = p + ggplot2::geom_vline(xintercept = vline)
        idx = findInterval(vline, x$t)
        if (idx >= 1 && idx <= nrow(x) && length(visible) > 0) {
            vline_data = data.frame(
                t = vline,
                name = factor(visible, levels = all_names),
                value = as.numeric(x[idx, visible])
            )
            p = p + ggplot2::geom_point(data = vline_data,
                ggplot2::aes(x = t, y = value, colour = name),
                shape = 4, size = 3, stroke = 0.5, show.legend = FALSE)
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
                ggplot2::aes(x = .data[[time_col]], y = value, colour = name), size = 2.5) +
            ggplot2::geom_point(data = overlay_long,
                ggplot2::aes(x = .data[[time_col]], y = value),
                shape = 1, colour = "black", size = 2.5, stroke = 0.4)
    }

    # Apply palette
    palette_fn = resolve_palette(palette)
    n_colours = length(all_names)
    p = p + ggplot2::scale_colour_manual(values = palette_fn(n_colours))
    if (use_linetypes) {
        p = p + ggplot2::scale_linetype_manual(values = grey_linetypes(n_colours))
    }

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
