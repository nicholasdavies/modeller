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
#' @param hide Optional character vector of compartment names to omit from the
#'   plot and from the initial conditions controls. Hidden compartments are
#'   still simulated using their default initial values.
#' @param max_display_rows Maximum rows shown in the Table tab (default
#'   5000). Full data is always available via CSV export.
#'
#' @return A Shiny app object (launched when printed).
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
#' if (interactive()) show_model(m)
#'
#' @export
show_model = function(model, data = NULL, hide = NULL, max_display_rows = 5000)
{
    init = model$init
    params = model$params
    extra_elements = model$shiny$ui
    hide = as.character(hide)

    # The interactive UI doesn't yet support vector-valued compartments or
    # vector-valued parameters. Bail out early with a clear message rather
    # than failing inside an inshiny widget.
    if (any(vapply(init, length, integer(1)) > 1L) ||
        any(vapply(params[setdiff(names(params), "time")], length, integer(1)) > 1L)) {
        stop("show_model() does not yet support models with vector-valued ",
            "compartments or parameters. Use run_model() and plot() directly.",
            call. = FALSE)
    }

    # Create Shiny app

    init_elements = list(shiny::tags$p(shiny::tags$strong("Initial conditions")))
    init_defaults = list()
    for (n in names(init)) {
        if (startsWith(n, cumulative_prefix)) next
        if (n %in% hide) next
        id = paste0("model_init_", n)
        init_defaults[[id]] = init[[n]]
        init_elements[[length(init_elements) + 1]] = inshiny::inline(n, "(0) = ",
            inshiny::inline_number(id, value = init[[n]],
                placeholder = init[[n]], min = 0, arrows = FALSE))
    }

    param_elements = list(shiny::tags$p(shiny::tags$strong("Parameters")))
    param_defaults = list()
    for (n in setdiff(names(params), "time")) {
        id = paste0("model_param_", n)
        param_defaults[[id]] = params[[n]]
        param_elements[[length(param_elements) + 1]] = inshiny::inline(n, " = ",
            inshiny::inline_number(id, value = params[[n]],
                placeholder = init[[n]], min = 0, arrows = FALSE))
    }

    # Reset button: clears any user changes to inputs above (init, params,
    # solver settings) and restores the values they had when show_model() was
    # called. Plot/UI options (palette, legend position, etc.) are unaffected.
    reset_defaults = c(init_defaults, param_defaults,
        model$shiny$defaults %||% list())
    reset_element = shiny::div(style = "margin-top: 15px;",
        shiny::actionButton("reset_inputs", "Reset to defaults",
            icon = shiny::icon("rotate-left"), class = "btn-sm"))

    # Colour palette choices
    palette_choices = c(
        "Okabe-Ito", "R4", "Tableau 10", "Dark2", "Set1", "Paired",
        "ggplot2", "viridis", "Grey"
    )

    # UI
    ui = bslib::page_sidebar(
        theme = bslib::bs_theme(version = 5, preset = "flatly", font_scale = 0.9),
        fillable = FALSE,
        shiny::tags$style(
            ".checkbox-inline { margin-right: 0.5em; }",
            ".alert-dismissible .btn-close { top: 0; transform: none; padding: 0.75rem; }"
        ),
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
                shiny::div(style = "display: flex; justify-content: flex-end; margin-bottom: 5px;",
                    inshiny::inline("Compare:  ",
                        inshiny::inline_switch("toggle_compare", value = FALSE,
                            on = "On", off = "Off",
                            meaning = "Save current model output to compare against future runs"))),
                shiny::plotOutput("main_plot", click = "main_plot_click"),
                shiny::uiOutput("recorded_plot_container"),
                shiny::uiOutput("message_box"),
                shiny::uiOutput("coords_box"),
                shiny::h5("Plot options", class = "mt-3"),
                shiny::uiOutput("compartment_toggles"),
                inshiny::inline("Plot title: ",
                    inshiny::inline_text("plot_title", value = NULL, placeholder = "(none)")),
                inshiny::inline("X label: ",
                    inshiny::inline_text("plot_xlab", value = "Time", placeholder = "(none)")),
                inshiny::inline("Y label: ",
                    inshiny::inline_text("plot_ylab", value = NULL, placeholder = "(none)")),
                inshiny::inline("Colour palette: ",
                    inshiny::inline_select("colour_palette", palette_choices, selected = "Dark2"),
                    " ", shiny::uiOutput("palette_preview", inline = TRUE)),
                inshiny::inline("Legend: ",
                    inshiny::inline_select("legend_position",
                        c("right", "left", "top", "bottom", "none"),
                        selected = "right")),
                shiny::br(),
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
            bslib::nav_panel("Load data",
                shiny::fileInput("data_upload", "Load CSV file", accept = ".csv"),
                shiny::uiOutput("data_time_col_ui"),
                shiny::tableOutput("data_preview")
            )
        ),
        sidebar = bslib::sidebar(title = "Model explorer", open = "always",
            init_elements, param_elements, extra_elements, reset_element)
    )

    server = function(input, output) {
        selected_time = shiny::reactiveVal()
        model_data = shiny::reactiveVal()
        imported_data = shiny::reactiveVal(data)
        last_plot = shiny::reactiveVal()
        model_messages = shiny::reactiveVal(NULL)
        compare_data = shiny::reactiveVal(NULL)

        # Compare toggle: snapshot current output when on, clear when off
        shiny::observeEvent(input$toggle_compare, {
            if (isTRUE(input$toggle_compare)) {
                compare_data(shiny::isolate(current_data()))
            } else {
                compare_data(NULL)
            }
        }, ignoreInit = TRUE)

        # Click handling: snap to nearest point in the data
        shiny::observeEvent(input$main_plot_click, {
            d = model_data()
            shiny::req(d)
            x = input$main_plot_click$x
            idx = which.min(abs(d$t - x))
            selected_time(d$t[idx])
        })

        shiny::observeEvent(input$recorded_plot_click, {
            d = model_data()
            shiny::req(d)
            x = input$recorded_plot_click$x
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

        shiny::observeEvent(input$clear_messages, {
            model_messages(NULL)
        })

        # Reset button: push the stored defaults back into each input.
        shiny::observeEvent(input$reset_inputs, {
            for (id in names(reset_defaults)) {
                inshiny::update_inline(id, value = reset_defaults[[id]])
            }
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
        shiny::outputOptions(output, "data_time_col_ui", suspendWhenHidden = FALSE)

        output$data_preview = shiny::renderTable({
            shiny::req(imported_data())
            head(imported_data(), 10)
        })

        # Cheap trigger reactive: depends on every input whose ID starts
        # with one of the model-related prefixes. Debouncing it defers the
        # expensive model run until inputs settle, so holding an arrow key
        # on a numeric input doesn't queue many runs. UI-only inputs (plot
        # title, palette, etc.) are excluded so changing them doesn't trigger
        # a model re-run.
        model_input_pattern = "^(model_init_|model_param_|ode_|ssa_|diff_)"
        model_inputs_trigger = shiny::reactive({
            for (n in grep(model_input_pattern, names(input), value = TRUE)) input[[n]]
            invisible(NULL)
        })
        debounced_inputs = shiny::debounce(model_inputs_trigger, 250)

        # Shared reactive: run model with current inputs (fires only after
        # the debounced trigger settles).
        current_data = shiny::reactive({
            debounced_inputs()
            shiny::isolate({
                init_list = lapply(names(init), function(n) {
                    if (startsWith(n, cumulative_prefix) || n %in% hide) init[[n]] else input[[paste0("model_init_", n)]]
                })
                names(init_list) = names(init)

                params_list = lapply(names(params), function(n) input[[paste0("model_param_", n)]])
                names(params_list) = names(params)

                model_messages(NULL)
                msgs = character(0)
                result = withCallingHandlers(
                    tryCatch(
                        model$shiny$run(model, init_list, params_list, input),
                        error = function(e) {
                            model_messages(list(type = "error", text = conditionMessage(e)))
                            NULL
                        }
                    ),
                    warning = function(w) {
                        msgs[[length(msgs) + 1L]] <<- conditionMessage(w)
                        model_messages(list(type = "warning", text = paste(msgs, collapse = "\n")))
                        invokeRestart("muffleWarning")
                    }
                )

                result
            })
        })

        # Identify columns produced by record() — present in data but not in
        # the model's compartments or incidence (renamed total_) columns.
        identify_recorded = function(d) {
            init_names = names(model$init)
            total_names = grep("^total_", init_names, value = TRUE)
            expected = c("t", setdiff(init_names, total_names),
                sub("^total_", "", total_names))
            setdiff(names(d), expected)
        }

        # Render compartment toggle checkboxes
        output$compartment_toggles = shiny::renderUI({
            d = current_data()
            shiny::req(d)
            series = setdiff(names(d), c("t", hide))
            prev = shiny::isolate(input$visible_series)
            selected = if (is.null(prev)) series else intersect(prev, series)
            shiny::checkboxGroupInput("visible_series", NULL,
                choices = series, selected = selected, inline = TRUE)
        })

        # Build main and (optional) recorded ggplot objects together so they
        # can be aligned via align_margins(); each is then rendered into its
        # own plotOutput so clicks resolve to data coords correctly.
        plot_pair = shiny::reactive({
            d = current_data()
            shiny::req(d)
            shiny::req(input$visible_series)

            imported = imported_data()
            time_col = if (!is.null(imported)) input$data_time_col %||% "t" else "t"
            if (!is.null(imported)) shiny::req(input$data_time_col)

            blank_to_null = function(x) if (is.null(x) || !nzchar(x)) NULL else x
            recorded = identify_recorded(d)
            visible_recorded = intersect(input$visible_series, recorded)
            visible_main = setdiff(input$visible_series, recorded)

            p_main = plot(d,
                series = visible_main,
                overlay = imported,
                palette = input$colour_palette %||% "Dark2",
                legend = input$legend_position %||% "right",
                time_col = time_col,
                vline = selected_time(),
                title = blank_to_null(input$plot_title),
                xlab = blank_to_null(input$plot_xlab),
                ylab = blank_to_null(input$plot_ylab),
                compare = compare_data())

            if (length(visible_recorded) == 0) {
                return(list(main = p_main, rec = NULL))
            }

            d_rec = d[, c("t", recorded), drop = FALSE]
            class(d_rec) = class(d)
            attr(d_rec, "dt") = attr(d, "dt")
            attr(d_rec, "geom") = attr(d, "geom")

            cmp = compare_data()
            cmp_rec = if (!is.null(cmp)) {
                cmp_cols = intersect(recorded, names(cmp))
                if (length(cmp_cols) > 0) {
                    x = cmp[, c("t", cmp_cols), drop = FALSE]
                    class(x) = class(cmp)
                    attr(x, "dt") = attr(cmp, "dt")
                    attr(x, "geom") = attr(cmp, "geom")
                    x
                } else NULL
            } else NULL

            # Continue palette after the colours the main plot uses
            n_main = length(setdiff(names(d), c("t", recorded)))
            n_total = n_main + length(recorded)
            base_fn = resolve_palette(input$colour_palette %||% "Dark2")
            offset_palette = function(n) base_fn(n_total)[(n_main + 1):(n_main + n)]

            p_rec = plot(d_rec,
                series = visible_recorded,
                palette = offset_palette,
                legend = input$legend_position %||% "right",
                vline = selected_time(),
                xlab = blank_to_null(input$plot_xlab),
                ylab = paste(visible_recorded, collapse = ", "),
                compare = cmp_rec)

            aligned = align_margins(list(p_main, p_rec), x = TRUE, y = FALSE)
            list(main = aligned[[1]], rec = aligned[[2]])
        })

        output$main_plot = shiny::renderPlot({
            pp = plot_pair()
            model_data(current_data())
            last_plot(pp$main)
            pp$main
        })

        # Conditional second plot for recorded quantities — only rendered
        # when at least one recorded series is selected.
        output$recorded_plot_container = shiny::renderUI({
            pp = plot_pair()
            if (is.null(pp$rec)) return(NULL)
            shiny::plotOutput("recorded_plot",
                click = "recorded_plot_click", height = "200px")
        })

        output$recorded_plot = shiny::renderPlot({
            pp = plot_pair()
            shiny::req(pp$rec)
            pp$rec
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

        # Message box: shows warnings/errors from model run
        output$message_box = shiny::renderUI({
            msg = model_messages()
            shiny::req(msg)
            alert_class = if (msg$type == "error") "alert-danger" else "alert-warning"
            prefix = if (msg$type == "error") "Error: " else "Warning: "
            shiny::div(class = paste("alert", alert_class, "alert-dismissible mb-2 mt-2 py-2"),
                style = "font-size: 0.9em; max-height: 200px; overflow-y: auto;",
                shiny::tags$span(shiny::tags$strong(prefix), msg$text),
                shiny::tags$button(type = "button", class = "btn-close",
                    onclick = "Shiny.setInputValue('clear_messages', Math.random(), {priority: 'event'});")
            )
        })


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
                style = "font-size: 0.9em; max-height: 200px; overflow-y: auto;",
                shiny::tags$span(style = "font-family: monospace; white-space: pre-wrap;", coord_text),
                shiny::tags$button(type = "button", class = "btn-close",
                    onclick = "Shiny.setInputValue('clear_coords', Math.random(), {priority: 'event'});")
            )
        })
    }

    app = shiny::shinyApp(ui, server)
    if (interactive()) print(app)
    invisible(app)
}
