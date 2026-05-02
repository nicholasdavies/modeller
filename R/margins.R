# Panel rects for ggplot2 plots
panel_rects_ggplot = function(plot)
{
    devsize = graphics::par("fin") * 25.4 # Device size in mm
    built = ggplot2::ggplot_build(plot)
    gtable = ggplot2::ggplot_gtable(built)
    to_mm = function(x) grid::convertUnit(x, "mm", valueOnly = TRUE)

    # With a set of lengths (widths/heights) from a gtable and a total length
    # (plot width/height in mm), replace null lengths with mm lengths
    expand_length = function(lengths, total_mm)
    {
        leeway_mm = total_mm - sum(to_mm(lengths))
        null_length = sum(sapply(unclass(lengths),
            function(x) if (x[[3]] == 5L) x[[1]] else 0))

        unit_class = class(lengths)
        lengths = unclass(lengths)
        for (l in seq_along(lengths)) {
            if (lengths[[l]][[3]] == 5L) { # 5L is code for null units
                lengths[[l]][[1]] = leeway_mm * lengths[[l]][[1]] / null_length
                lengths[[l]][[3]] = 7L # change to mm
            }
        }
        class(lengths) = unit_class

        return (lengths)
    }

    # Get width and height in mm for gtable as realized plot
    w = to_mm(expand_length(lengths = gtable$widths, total_mm = devsize[1]))
    h = to_mm(expand_length(lengths = gtable$heights, total_mm = devsize[2]))

    # Create data.frame to return
    layout = as.data.frame(gtable[[2]])
    p = layout[grepl("^panel", layout$name), ]
    p$x = sapply(1:nrow(p), function(i) sum(w[seq_len(p$l[i] - 1)]) / devsize[1])
    p$y = sapply(1:nrow(p), function(i) 1 - sum(h[seq_len(p$t[i])]) / devsize[2])
    p$w = sapply(1:nrow(p), function(i) sum(w[p$l[i]:p$r[i]]) / devsize[1])
    p$h = sapply(1:nrow(p), function(i) sum(h[p$t[i]:p$b[i]]) / devsize[2])

    # Add coordinate information
    # This is more complicated than seems necessary because ggplot doesn't
    # necessarily seem to name the panels correctly, e.g. panel-X-Y may not be
    # either the panel at row X column Y, nor the panel at row Y column X.
    panel_grobids = which(grepl("^panel", gtable[[2]]$name)) # grob indices of panels
    panel_ggnames = gtable[[2]]$name[panel_grobids] # ggplot2 names of panels
    panel_grobnames = sapply(panel_grobids, function(i) gtable[[1]][[i]]$name) # names of grobs for each panel
    panel_pids = as.integer(sapply(regmatches(panel_grobnames,
        regexec("^panel-([0-9]+)\\.", panel_grobnames)), `[`, 2)) # panel ids (1 to N) for each panel

    # Add panel id to each row of the p data frame (NA for not a valid panel...)
    p$panel_id = panel_pids[match(p$name, panel_ggnames)]
    p = p[!is.na(p$panel_id), ]

    # Add coordinate limits:
    p$xmin = sapply(p$panel_id, function(i) built$layout$panel_params[[i]]$x.range[1])
    p$xmax = sapply(p$panel_id, function(i) built$layout$panel_params[[i]]$x.range[2])
    p$ymin = sapply(p$panel_id, function(i) built$layout$panel_params[[i]]$y.range[1])
    p$ymax = sapply(p$panel_id, function(i) built$layout$panel_params[[i]]$y.range[2])

    # Add layout information
    p_id_order = match(p$panel_id, built$layout$layout$PANEL)
    p$row = built$layout$layout[p_id_order, "ROW"]
    p$col = built$layout$layout[p_id_order, "COL"]

    # Add label
    lab0 = which(names(built$layout$layout) == "COL")
    lab1 = which(names(built$layout$layout) == "SCALE_X")
    if (lab1 == lab0 + 1) {
        p$label = NA_character_
    } else if (lab1 == lab0 + 2) {
        p$label = as.character(built$layout$layout[p_id_order, lab0 + 1])
    } else {
        p$label = paste(
            built$layout$layout[p_id_order, lab0 + 1],
            built$layout$layout[p_id_order, lab0 + 2],
            sep = "~")
    }

    # Tidy data frame
    p = p[, c("panel_id", "row", "col", "label", "x", "y", "w", "h", "xmin", "xmax", "ymin", "ymax", "name")]
    names(p)[13] = "grobname"
    p = p[order(p$panel_id), ]
    rownames(p) = 1:nrow(p)

    return (p)
}

# Convert x to normalised plot coordinates
to_npc = function(x)
{
    grid::convertUnit(x, "npc", valueOnly = TRUE)
}

# With a set of lengths (widths/heights) from a gtable, return lengths
# in normalised parent coordinates of the plot margins
get_margins = function(lengths)
{
    nulls = sapply(unclass(lengths), function(x) x[[3]] == 5L)
    if (!any(nulls)) {
        stop("Could not retrieve plot margins.")
    }
    first = which(nulls)[1]
    last = which(nulls)[sum(nulls)]
    length0 = sum(to_npc(lengths[1:(first - 1)]))
    length1 = sum(to_npc(lengths[(last + 1):length(lengths)]))
    return (c(length0, length1))
}

# Align margins of a set of ggplot2 plots
#
# Computes the maximum non-panel margin (in mm) across the input plots, then
# pads each plot's theme(plot.margin) by the difference so panels line up.
# Returns ggplot objects (not ggdraws), so the original coordmap is preserved
# and downstream tools such as shiny click handlers see data coordinates.
#
# gl  a list of ggplot2 plots.
# x   if TRUE, align x axes
# y   if TRUE, align y axes
#
# returns a list of aligned plots
align_margins = function(gl, x, y)
{
    devsize_mm = graphics::par("fin") * 25.4
    to_mm = function(u) grid::convertUnit(u, "mm", valueOnly = TRUE)

    # Expand a unit vector (which may include "null" units) to mm given a
    # known total length in mm.
    expand_lengths_mm = function(lengths, total_mm) {
        n = length(lengths)
        nulls = grid::unitType(lengths) == "null"
        if (!any(nulls)) return (to_mm(lengths))
        mm = numeric(n)
        null_weights = numeric(n)
        for (i in seq_len(n)) {
            if (nulls[i]) {
                null_weights[i] = as.numeric(lengths[i])
            } else {
                mm[i] = to_mm(lengths[i])
            }
        }
        null_total = sum(null_weights)
        leeway_mm = total_mm - sum(mm)
        for (i in seq_len(n)) {
            if (nulls[i]) mm[i] = leeway_mm * null_weights[i] / null_total
        }
        mm
    }

    # Build grobs and find first panel's row/col + non-panel margins in mm
    grobs = lapply(gl, cowplot::as_grob)
    panel_pos = lapply(grobs, function(g) {
        layout = as.data.frame(g[[2]])
        p = layout[grepl("^panel", layout$name), ][1, ]
        c(l = p$l, r = p$r, t = p$t, b = p$b)
    })
    h_margins = mapply(function(g, pp) {
        widths_mm = expand_lengths_mm(g$widths, devsize_mm[1])
        c(left  = sum(widths_mm[seq_len(pp["l"] - 1)]),
          right = sum(widths_mm[(pp["r"] + 1):length(widths_mm)]))
    }, grobs, panel_pos, SIMPLIFY = FALSE)
    v_margins = mapply(function(g, pp) {
        heights_mm = expand_lengths_mm(g$heights, devsize_mm[2])
        c(top    = sum(heights_mm[seq_len(pp["t"] - 1)]),
          bottom = sum(heights_mm[(pp["b"] + 1):length(heights_mm)]))
    }, grobs, panel_pos, SIMPLIFY = FALSE)

    max_l = max(sapply(h_margins, `[`, "left"))
    max_r = max(sapply(h_margins, `[`, "right"))
    max_t = max(sapply(v_margins, `[`, "top"))
    max_b = max(sapply(v_margins, `[`, "bottom"))

    default_pm = ggplot2::theme_get()$plot.margin
    if (is.null(default_pm)) default_pm = ggplot2::margin()

    lapply(seq_along(gl), function(i) {
        plot = gl[[i]]
        pm = plot$theme$plot.margin
        if (is.null(pm)) pm = default_pm
        pm_mm = to_mm(pm) # t, r, b, l

        new_l = pm_mm[4]; new_r = pm_mm[2]
        new_t = pm_mm[1]; new_b = pm_mm[3]
        if (x) {
            new_l = new_l + (max_l - h_margins[[i]]["left"])
            new_r = new_r + (max_r - h_margins[[i]]["right"])
        }
        if (y) {
            new_t = new_t + (max_t - v_margins[[i]]["top"])
            new_b = new_b + (max_b - v_margins[[i]]["bottom"])
        }

        plot + ggplot2::theme(plot.margin = ggplot2::margin(
            t = new_t, r = new_r, b = new_b, l = new_l, unit = "mm"))
    })
}
