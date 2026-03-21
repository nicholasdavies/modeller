#' @export
ode_model = function(init, params, model, times, ...)
{
    # Put model function into correct form
    # TODO check function's signature
    formals(model) = alist(t =, .state =, .params =)
    body(model) = rlang::expr({
        .x <- with(as.list(c(.state, .params)), !!body(model))
        return (list(unlist(.x)))
    })

    # Process times
    times$start = times$start %||% 0
    times$stop = times$stop %||% 100
    times$duration = times$duration %||% (times$stop - times$start)
    times$step = times$step %||% 1

    # Extra elements for model configuration
    ode_methods = c("lsoda", "lsode", "lsodes", "lsodar", "vode", "daspk",
           "euler", "rk4", "ode23", "ode45", "radau",
           "bdf", "bdf_d", "adams", "impAdams", "impAdams_d")
    extra_elements = list(
        shiny::h5("ODE settings"),
        inshiny::inline("Method: ", inshiny::inline_select("ode_method",
            choices = ode_methods, meaning = "ODE method")),
        inshiny::inline("Start time: ",
            inshiny::inline_number("ode_t0", value = times$start, min = 0,
                placeholder = times$start, arrows = FALSE, meaning = "Start time")),
        inshiny::inline("Duration: ",
            inshiny::inline_number("ode_tt", value = times$duration, min = 0,
                placeholder = times$duration, arrows = FALSE, meaning = "Duration")),
        inshiny::inline("Time step: ",
            inshiny::inline_number("ode_dt", value = times$step, min = 0.01, step = 0.01,
                placeholder = times$step, arrows = FALSE, meaning = "Time step"))

    )

    # Solver: returns data frame with t + compartment columns
    run_model = function(init_list, params_list, times, input) {
        tval = seq(input$ode_t0, input$ode_t0 + input$ode_tt, input$ode_dt)
        sol = deSolve::ode(y = unlist(init_list), times = tval,
            func = model, parms = params_list, method = input$ode_method)
        data = as.data.frame(sol)
        names(data)[1] = "t"
        attr(data, "dt") = input$ode_dt
        data
    }

    model_app(run_model, init, params, times, extra_elements, ...)
}
