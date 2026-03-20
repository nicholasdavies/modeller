model_app = function(run_model, init, params, times)
{
    # Process times
    time_start = times$start %||% 0
    time_stop = times$stop %||% 100
    time_step = times$step %||% 1
    times = seq(time_start, time_stop, time_step)

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

    # UI
    ui = bslib::page_sidebar(
        theme = bslib::bs_theme(version = 5, preset = "flatly", font_scale = 0.9),
        fillable = FALSE,
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
        shiny::plotOutput("main_plot", click = "main_plot_click"),
        shiny::uiOutput("coords_box"),
        sidebar = bslib::sidebar(title = "Settings", open = "always",
            init_elements, shiny::br(), param_elements)
    )

    server = function(input, output) {
        # Click handling: snap to nearest time step, show compartment values
        snap_to_time = function(x) {
            idx = which.min(abs(times - x))
            times[idx]
        }

        selected_time = shiny::reactiveVal()
        model_data = shiny::reactiveVal()

        shiny::observeEvent(input$main_plot_click, {
            x = input$main_plot_click$x
            if (x >= time_start && x <= time_stop) {
                selected_time(snap_to_time(x))
            }
        })

        # Arrow keys step through time when no input field has focus
        shiny::observeEvent(input$arrow_key, {
            current = selected_time()
            shiny::req(current)
            if (input$arrow_key$key == "ArrowRight") {
                new_t = min(current + time_step, time_stop)
            } else {
                new_t = max(current - time_step, time_start)
            }
            selected_time(new_t)
        })

        shiny::observeEvent(input$clear_coords, {
            selected_time(NULL)
        })

        output$main_plot = shiny::renderPlot({
            # Get initial conditions and parameters
            init_list = lapply(names(init), function(n) input[[paste0("model_init_", n)]])
            names(init_list) = names(init)

            params_list = lapply(names(params), function(n) input[[paste0("model_param_", n)]])
            names(params_list) = names(params)

            # Run model
            data = run_model(init_list, params_list, times)
            model_data(data)

            # Plot results
            plot_data = tidyr::pivot_longer(data, cols = -1)
            plot_data$name = factor(plot_data$name, unique(plot_data$name))

            ggplot2::ggplot(plot_data) +
                ggplot2::geom_line(ggplot2::aes(x = t, y = value, colour = name)) +
                ggplot2::geom_vline(xintercept = selected_time())
        })

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
