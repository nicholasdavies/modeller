#' Launch interactive model explorer
#'
#' Opens a Shiny app for interactively exploring a model. Provides controls
#' for adjusting initial conditions, parameters, and solver settings, with
#' live-updating plots and data tables. Observed data can be overlaid for
#' visual comparison.
#'
#' @param model A model object created by [ode_model()] or [ssa_model()].
#' @param data Optional data frame to overlay on the plot. Should contain a
#'   time column and one or more columns matching compartment or incidence
#'   names. Can also be loaded interactively via the Data tab.
#' @param max_display_rows Maximum rows shown in the Table tab (default
#'   5000). Full data is always available via CSV export.
#'
#' @return A Shiny app object (launched when printed).
#'
#' @examples
#' m = ode_model(
#'     init = list(S = 999, I = 1, R = 0),
#'     params = list(beta = 0.3, gamma = 0.1),
#'     model = function() list(-beta * S * I, beta * S * I - gamma * I, gamma * I),
#'     times = list(duration = 50)
#' )
#' if (interactive()) show_model(m)
#'
#' @export
show_model = function(model, data = NULL, max_display_rows = 5000)
{
    init = model$init
    params = model$params
    extra_elements = model$shiny$ui

    # Create Shiny app

    init_elements = list(shiny::tags$p(shiny::tags$strong("Initial conditions")))
    for (n in names(init)) {
        if (startsWith(n, "total_")) next
        id = paste0("model_init_", n)
        init_elements[[length(init_elements) + 1]] = inshiny::inline(n, "(0) = ",
            inshiny::inline_number(id, value = init[[n]],
                placeholder = init[[n]], min = 0, arrows = FALSE))
    }

    param_elements = list(shiny::tags$p(shiny::tags$strong("Parameters")))
    for (n in names(params)) {
        id = paste0("model_param_", n)
        param_elements[[length(param_elements) + 1]] = inshiny::inline(n, " = ",
            inshiny::inline_number(id, value = params[[n]],
                placeholder = init[[n]], min = 0, arrows = FALSE))
    }

    # Colour palette choices
    palette_choices = c(
        "Okabe-Ito", "R4", "Tableau 10", "Dark2", "Set1", "Paired",
        "ggplot2", "viridis", "magma", "plasma", "inferno"
    )

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
            bslib::nav_panel("Plot",
                shiny::plotOutput("main_plot", click = "main_plot_click"),
                shiny::uiOutput("coords_box"),
                shiny::h5("Plot options", class = "mt-3"),
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
                shiny::div(shiny::downloadButton("export_plot", "Export plot",
                    class = "btn-sm"))
            ),
            bslib::nav_panel("Table",
                shiny::div(style = "margin: 10px 0px 10px 0px",
                    shiny::downloadButton("export_csv", "Export CSV")),
                shiny::tableOutput("model_table"),
                shiny::uiOutput("table_cap")
            ),
            bslib::nav_panel("Data",
                shiny::fileInput("data_upload", "Load CSV file", accept = ".csv"),
                shiny::uiOutput("data_time_col_ui"),
                shiny::tableOutput("data_preview")
            )
        ),
        sidebar = bslib::sidebar(title = "Model explorer", open = "always",
            init_elements, param_elements, extra_elements)
    )

    server = function(input, output) {
        selected_time = shiny::reactiveVal()
        model_data = shiny::reactiveVal()
        imported_data = shiny::reactiveVal(data)
        last_plot = shiny::reactiveVal()

        # Click handling: snap to nearest point in the data
        shiny::observeEvent(input$main_plot_click, {
            d = model_data()
            shiny::req(d)
            x = input$main_plot_click$x
            idx = which.min(abs(d$t - x))
            selected_time(d$t[idx])
        })

        # Arrow keys step to next/previous data point
        shiny::observeEvent(input$arrow_key, {
            current = selected_time()
            shiny::req(current)
            d = model_data()
            shiny::req(d)
            idx = which.min(abs(d$t - current))
            if (input$arrow_key$key == "ArrowRight") {
                idx = min(idx + 1L, nrow(d))
            } else {
                idx = max(idx - 1L, 1L)
            }
            selected_time(d$t[idx])
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
            init_list = lapply(names(init), function(n) {
                if (startsWith(n, "total_")) init[[n]] else input[[paste0("model_init_", n)]]
            })
            names(init_list) = names(init)

            params_list = lapply(names(params), function(n) input[[paste0("model_param_", n)]])
            names(params_list) = names(params)

            model$shiny$run(model, init_list, params_list, input)
        })

        # Render compartment toggle checkboxes
        output$compartment_toggles = shiny::renderUI({
            d = current_data()
            shiny::req(d)
            total_cols = grep("^total_", names(d), value = TRUE)
            series = setdiff(names(d), c("t", total_cols))
            inc_names = sub("^total_", "", total_cols)
            series = c(series, inc_names)

            shiny::checkboxGroupInput("visible_series", NULL,
                choices = series, selected = series, inline = TRUE)
        })

        output$main_plot = shiny::renderPlot({
            d = current_data()
            model_data(d)

            imported = imported_data()
            time_col = if (!is.null(imported)) input$data_time_col %||% "t" else "t"
            if (!is.null(imported)) shiny::req(input$data_time_col)

            shiny::req(input$visible_series)

            p = plot(d,
                series = input$visible_series,
                overlay = imported,
                palette = input$colour_palette %||% "Set1",
                legend = input$legend_position %||% "right",
                time_col = time_col,
                vline = selected_time())
            last_plot(p)
            p
        })

        output$palette_preview = shiny::renderUI({
            palette = input$colour_palette
            shiny::req(palette)
            n = 6
            palette_fn = resolve_palette(palette)
            cols = palette_fn(n)
            swatches = lapply(cols, function(col) {
                shiny::tags$span(style = paste0(
                    "display:inline-block;width:12px;height:12px;",
                    "border-radius:2px;background:", col, ";margin-right:2px;",
                    "vertical-align:middle;"))
            })
            shiny::tagList(swatches)
        })

        output$model_table = shiny::renderTable({
            d = current_data()
            if (nrow(d) > max_display_rows) d = d[seq_len(max_display_rows), ]
            # Format each column: 0 decimal places if integer-valued, 2 otherwise
            for (col in names(d)) {
                if (all(d[[col]] == round(d[[col]]), na.rm = TRUE)) {
                    d[[col]] = formatC(d[[col]], format = "f", digits = 0, big.mark = ",")
                } else {
                    d[[col]] = formatC(d[[col]], format = "f", digits = 2, big.mark = ",")
                }
            }
            d
        })

        output$table_cap = shiny::renderUI({
            d = current_data()
            if (nrow(d) > max_display_rows) {
                shiny::p(style = "color: #888; font-style: italic;",
                    paste0("Display capped at ", format(max_display_rows, big.mark = ","),
                        " rows (", format(nrow(d), big.mark = ","), " total). ",
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

        # Plot export
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
            d = model_data()
            shiny::req(d)

            row = d[d$t == t_sel, , drop = FALSE]
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
