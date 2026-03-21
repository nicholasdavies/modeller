model_app = function(run_model, init, params, times, extra_elements = NULL,
    max_display_rows = 5000)
{
    # Create Shiny app

    # TODO make sure init is a list . . .
    init_elements = list(shiny::h5("Initial conditions"))
    for (n in names(init)) {
        id = paste0("model_init_", n)
        init_elements[[length(init_elements) + 1]] = inshiny::inline(n, "(0) = ",
            inshiny::inline_number(id, value = init[[n]],
                placeholder = init[[n]], min = 0, arrows = FALSE))
    }

    # TODO make sure params is a list . . .
    param_elements = list(shiny::h5("Parameters"))
    for (n in names(params)) {
        id = paste0("model_param_", n)
        param_elements[[length(param_elements) + 1]] = inshiny::inline(n, " = ",
            inshiny::inline_number(id, value = params[[n]],
                placeholder = init[[n]], min = 0, arrows = FALSE))
    }

    # Colour palette choices
    palette_choices = c(
        "Okabe-Ito", "R4", "Tableau 10", "Dark2", "Set1", "Paired",
        "viridis", "magma", "plasma", "inferno"
    )
    viridis_palettes = c("viridis", "magma", "plasma", "inferno")
    base_palettes = c("Okabe-Ito", "R4", "Tableau 10")

    # UI
    ui = bslib::page_sidebar(
        theme = bslib::bs_theme(version = 5, preset = "flatly", font_scale = 0.9),
        fillable = FALSE,
        shiny::tags$style(".checkbox-inline { margin-right: 0.5em; }"),
        # Arrow key listener (only when no input field is focused)
        shiny::tags$script(shiny::HTML("
            $(document).on('keydown', function(e) {
                if (e.key === 'ArrowLeft' || e.key === 'ArrowRight') {
                    var el = document.activeElement;
                    if (!el || el.tagName === 'BODY' || el.id === 'main_plot' || el.closest('#main_plot')) {
                        Shiny.setInputValue('arrow_key', {key: e.key, rand: Math.random()}, {priority: 'event'});
                        e.preventDefault();
                    }
                }
            });
        ")),
        bslib::navset_underline(
            bslib::nav_panel("Graph",
                shiny::plotOutput("main_plot", click = "main_plot_click"),
                shiny::uiOutput("coords_box"),
                bslib::card(class = "mt-2",
                    bslib::card_header("Graph options"),
                    bslib::card_body(
                        shiny::uiOutput("compartment_toggles"),
                        inshiny::inline("Colour palette: ",
                            inshiny::inline_select("colour_palette", palette_choices, selected = "Set1"),
                            " ", shiny::uiOutput("palette_preview", inline = TRUE)),
                        inshiny::inline("Legend: ",
                            inshiny::inline_select("legend_position",
                                c("right", "left", "top", "bottom", "none"),
                                selected = "right")),
                        inshiny::inline("Export as ",
                            inshiny::inline_select("export_format", c("PDF", "PNG"), selected = "PDF"),
                            ", ",
                            inshiny::inline_number("export_width", value = 20, min = 1, arrows = FALSE),
                            " \u00d7 ",
                            inshiny::inline_number("export_height", value = 12, min = 1, arrows = FALSE),
                            " cm"),
                        shiny::div(shiny::downloadButton("export_plot", "Export graph",
                            class = "btn-sm"))
                    )
                )
            ),
            bslib::nav_panel("Table",
                shiny::div(style = "margin: 10px 0px 10px 0px",
                    shiny::downloadButton("export_csv", "Export CSV")),
                shiny::tableOutput("model_table"),
                shiny::uiOutput("table_cap")
            ),
            bslib::nav_panel("Data",
                shiny::fileInput("data_upload", "Upload CSV file", accept = ".csv"),
                shiny::uiOutput("data_time_col_ui"),
                shiny::tableOutput("data_preview")
            )
        ),
        sidebar = bslib::sidebar(title = "Settings", open = "always",
            init_elements, shiny::br(), param_elements, shiny::br(),
            extra_elements)
    )

    server = function(input, output) {
        selected_time = shiny::reactiveVal()
        model_data = shiny::reactiveVal()
        imported_data = shiny::reactiveVal()
        last_plot = shiny::reactiveVal()

        # Click handling: snap to nearest point in the data
        shiny::observeEvent(input$main_plot_click, {
            data = model_data()
            shiny::req(data)
            x = input$main_plot_click$x
            idx = which.min(abs(data$t - x))
            selected_time(data$t[idx])
        })

        # Arrow keys step to next/previous data point
        shiny::observeEvent(input$arrow_key, {
            current = selected_time()
            shiny::req(current)
            data = model_data()
            shiny::req(data)
            idx = which.min(abs(data$t - current))
            if (input$arrow_key$key == "ArrowRight") {
                idx = min(idx + 1L, nrow(data))
            } else {
                idx = max(idx - 1L, 1L)
            }
            selected_time(data$t[idx])
        })

        shiny::observeEvent(input$clear_coords, {
            selected_time(NULL)
        })

        # Data import
        shiny::observeEvent(input$data_upload, {
            imported_data(data.table::fread(input$data_upload$datapath, data.table = FALSE))
        })

        output$data_time_col_ui = shiny::renderUI({
            shiny::req(imported_data())
            cols = names(imported_data())
            # Default to column named "t" or "time" (case-insensitive)
            default = cols[tolower(cols) %in% c("t", "time")]
            if (length(default) == 0) default = cols[1]
            shiny::selectInput("data_time_col", "Time column", choices = cols, selected = default[1])
        })

        output$data_preview = shiny::renderTable({
            shiny::req(imported_data())
            head(imported_data(), 10)
        })

        # Shared reactive: run model with current inputs
        current_data = shiny::reactive({
            init_list = lapply(names(init), function(n) input[[paste0("model_init_", n)]])
            names(init_list) = names(init)

            params_list = lapply(names(params), function(n) input[[paste0("model_param_", n)]])
            names(params_list) = names(params)

            run_model(init_list, params_list, times, input)
        })

        # All series names (compartments + incidence), updated when model changes
        all_series = shiny::reactiveVal()

        # Render compartment toggle checkboxes
        output$compartment_toggles = shiny::renderUI({
            data = current_data()
            shiny::req(data)
            total_cols = grep("^total_", names(data), value = TRUE)
            series = setdiff(names(data), c("t", total_cols))
            inc_names = sub("^total_", "", total_cols)
            series = c(series, inc_names)
            all_series(series)

            shiny::checkboxGroupInput("visible_series", NULL,
                choices = series, selected = series, inline = TRUE)
        })

        output$main_plot = shiny::renderPlot({
            data = current_data()
            model_data(data)

            # Separate total_ columns for incidence overlay
            total_cols = grep("^total_", names(data), value = TRUE)
            plot_cols = setdiff(names(data), total_cols)
            main_data = data[plot_cols]

            # Determine visible series
            visible = input$visible_series %||% character(0)

            # Plot results
            plot_data = tidyr::pivot_longer(main_data, cols = -1)
            plot_data$name = factor(plot_data$name, unique(plot_data$name))
            plot_data = plot_data[plot_data$name %in% visible, ]

            compartment_cols = setdiff(names(main_data), "t")
            all_integer = all(main_data[compartment_cols] == round(main_data[compartment_cols]))
            geom = if (all_integer) ggplot2::geom_step else ggplot2::geom_line

            p = ggplot2::ggplot(plot_data) +
                geom(ggplot2::aes(x = t, y = value, colour = name), linewidth = 1) +
                ggplot2::geom_vline(xintercept = selected_time()) +
                ggplot2::labs(x = "Time", y = NULL, colour = NULL) +
                nice_scales +
                cowplot::theme_cowplot(font_size = 12) +
                ggplot2::theme(
                    axis.line = ggplot2::element_line(linewidth = 0.3),
                    axis.ticks = ggplot2::element_line(linewidth = 0.3),
                    legend.position = input$legend_position %||% "right"
                )

            # Overlay incidence from total_ columns
            if (length(total_cols) > 0) {
                mat = as.matrix(data[c("t", total_cols)])
                grid_step = if (attr(data, "dt") <= 1) 1 else attr(data, "dt")
                grid = seq(min(mat[, 1]), max(mat[, 1]), by = grid_step)
                summ = grid_summarize(mat, grid)
                inc_data = data.frame(t = grid)
                for (col in total_cols) {
                    new_name = sub("^total_", "", col)
                    vals = summ[, col] / grid_step
                    inc_data[[new_name]] = c(vals[1], diff(vals))
                }
                inc_long = tidyr::pivot_longer(inc_data, cols = -1)
                inc_long = inc_long[inc_long$name %in% visible, ]
                all_levels = union(levels(plot_data$name), unique(inc_long$name))
                plot_data$name = factor(plot_data$name, levels = all_levels)
                inc_long$name = factor(inc_long$name, levels = all_levels)
                if (nrow(inc_long) > 0) {
                    p = p + geom(data = inc_long,
                        ggplot2::aes(x = t, y = value, colour = name), linewidth = 1)
                }
            }

            # Overlay imported data
            imported = imported_data()
            if (!is.null(imported)) {
                time_col = input$data_time_col
                shiny::req(time_col)
                other_cols = setdiff(names(imported), time_col)
                overlay = tidyr::pivot_longer(imported, cols = dplyr::all_of(other_cols),
                    names_to = "name", values_to = "value")
                overlay$name = factor(overlay$name, levels = levels(plot_data$name))
                overlay = overlay[overlay$name %in% visible, ]
                p = p + ggplot2::geom_point(data = overlay,
                    ggplot2::aes(x = .data[[time_col]], y = value, colour = name),
                    size = 2.5) +
                    ggplot2::geom_point(data = overlay,
                    ggplot2::aes(x = .data[[time_col]], y = value),
                    shape = 1, colour = "black", size = 2.5, stroke = 0.4)
            }

            # Apply colour palette
            palette = input$colour_palette
            if (!is.null(palette)) {
                if (palette %in% viridis_palettes) {
                    p = p + ggplot2::scale_colour_viridis_d(option = palette)
                } else if (palette %in% base_palettes) {
                    cols = unname(palette.colors(palette = palette))
                    p = p + ggplot2::scale_colour_manual(values = cols)
                } else {
                    p = p + ggplot2::scale_colour_brewer(palette = palette)
                }
            }

            last_plot(p)
            p
        })

        output$palette_preview = shiny::renderUI({
            palette = input$colour_palette
            shiny::req(palette)
            n = 6
            if (palette %in% viridis_palettes) {
                cols = scales::viridis_pal(option = palette)(n)
            } else if (palette %in% base_palettes) {
                cols = unname(palette.colors(n, palette = palette))
            } else {
                cols = RColorBrewer::brewer.pal(n, palette)
            }
            swatches = lapply(cols, function(col) {
                shiny::tags$span(style = paste0(
                    "display:inline-block;width:12px;height:12px;",
                    "border-radius:2px;background:", col, ";margin-right:2px;",
                    "vertical-align:middle;"))
            })
            shiny::tagList(swatches)
        })

        output$model_table = shiny::renderTable({
            data = current_data()
            if (nrow(data) > max_display_rows) data = data[seq_len(max_display_rows), ]
            # Format each column: 0 decimal places if integer-valued, 2 otherwise
            for (col in names(data)) {
                if (all(data[[col]] == round(data[[col]]), na.rm = TRUE)) {
                    data[[col]] = formatC(data[[col]], format = "f", digits = 0, big.mark = ",")
                } else {
                    data[[col]] = formatC(data[[col]], format = "f", digits = 2, big.mark = ",")
                }
            }
            data
        })

        output$table_cap = shiny::renderUI({
            data = current_data()
            if (nrow(data) > max_display_rows) {
                shiny::p(style = "color: #888; font-style: italic;",
                    paste0("Display capped at ", format(max_display_rows, big.mark = ","),
                        " rows (", format(nrow(data), big.mark = ","), " total). ",
                        "Use Export CSV for full data."))
            }
        })

        output$export_csv = shiny::downloadHandler(
            filename = function() {
                paste0("model_output_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
            },
            content = function(file) {
                data.table::fwrite(current_data(), file)
            }
        )

        # Graph export
        output$export_plot = shiny::downloadHandler(
            filename = function() {
                ext = tolower(input$export_format)
                paste0("model_plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".", ext)
            },
            content = function(file) {
                ggplot2::ggsave(file, plot = last_plot(),
                    width = input$export_width, height = input$export_height, units = "cm")
            }
        )

        # Coordinates box: shows compartment values at selected time
        output$coords_box = shiny::renderUI({
            t_sel = selected_time()
            shiny::req(t_sel)
            data = model_data()
            shiny::req(data)

            row = data[data$t == t_sel, , drop = FALSE]
            if (nrow(row) == 0) return(NULL)

            compartments = names(row)[names(row) != "t"]
            vals = paste0(compartments, " = ", signif(as.numeric(row[1, compartments]), 4))
            coord_text = paste0("t = ", t_sel, ":\n", paste(vals, collapse = ",\n"))

            shiny::div(class = "alert alert-secondary alert-dismissible mb-2 mt-2 py-2",
                style = "font-size: 0.9em;",
                shiny::tags$span(style = "font-family: monospace; white-space: pre-wrap;", coord_text),
                shiny::tags$button(type = "button", class = "btn-close",
                    onclick = "Shiny.setInputValue('clear_coords', Math.random(), {priority: 'event'});")
            )
        })
    }

    shiny::shinyApp(ui, server)
}
