#' Create a deterministic ODE model
#'
#' Defines a compartmental model solved with ordinary differential equations
#' via [deSolve::ode()]. At each time point, the equations function
#' calculates and returns the rate of change (derivative) of every
#' compartment. For example, if 10 individuals per unit time are leaving
#' compartment `S`, the equations should return `dS = -10`.
#'
#' The equations function is written as a plain R function with no arguments.
#' Inside the function body, compartment names (from `init`) and parameter
#' names (from `params`) can be used as ordinary variables — the package
#' makes them available automatically.
#'
#' @param init Named list of initial compartment values (e.g.
#'   `list(S = 999, I = 1, R = 0)`).
#' @param params Named list of parameter values (e.g.
#'   `list(beta = 0.3, gamma = 0.1)`).
#' @param equations A function with no arguments. Inside the function body,
#'   write equations using compartment and parameter names as if they were
#'   ordinary variables. Return a list giving the derivative (rate of change)
#'   of each compartment, in the same order as `init`. Compartments whose names
#'   start with `total_` are treated as cumulative counters; the corresponding
#'   incidence is computed automatically when plotting.
#' @param times Named list controlling the time span. Accepts `start`,
#'   `stop`, `duration`, and `step`. At minimum, provide `duration` or
#'   `stop`. Defaults: `start = 0`, `stop = 100`, `step = 1`.
#' @param options Named list of solver options. For ODE models the only
#'   option is `method`, the integration method passed to [deSolve::ode()].
#'   Default: `list(method = "lsoda")`. Can be overridden in [run_model()].
#'
#' @return An object of class `c("ode_model", "modeller")`. Use
#'   [run_model()] to simulate, [plot()][plot.model_result] to visualise,
#'   and [show_model()] to launch the interactive app.
#'
#' @examples
#' m = ode_model(
#'     init = list(S = 999, I = 1, R = 0),
#'     params = list(beta = 0.3, gamma = 0.1),
#'     equations = function() {
#'         N = S + I + R
#'         list(
#'             dS = -beta * I / N * S,
#'             dI =  beta * I / N * S - gamma * I,
#'             dR =  gamma * I
#'         )
#'     },
#'     times = list(duration = 100)
#' )
#' m
#' result = run_model(m)
#' plot(result)
#'
#' @export
ode_model = function(init, params, equations, times, options = list())
{
    validate_inputs(init, params, equations, times)

    # Put equations function into correct form
    formals(equations) = alist(t =, .state =, .params =)
    eq = body(equations)
    eq = elixir::expr_replace(eq, { d(`.A:name`) <- ..B },
        { .dlist[[as.character(quote(.A))]] <<- ..B })
    eq = elixir::expr_replace(eq, { d(`.A:name`) = ..B },
        { .dlist[[as.character(quote(.A))]] <<- ..B })

    body(equations) = rlang::expr({
        .dlist = init
        with(c(.state, .params), !!eq)
        return (list(unlist(.dlist)))
    })

    # Process times
    times = resolve_times(times, has_step = TRUE)

    # Shiny UI elements for solver settings
    ode_methods = c("lsoda", "lsode", "lsodes", "lsodar", "vode", "daspk",
           "euler", "rk4", "ode23", "ode45", "radau",
           "bdf", "bdf_d", "adams", "impAdams", "impAdams_d")
    default_options = modifyList(list(method = "lsoda"), options)
    shiny_ui = list(
        shiny::tags$p(shiny::tags$strong("ODE settings")),
        inshiny::inline("Method: ", inshiny::inline_select("ode_method",
            choices = ode_methods, selected = default_options$method,
            meaning = "ODE method")),
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

    # Shiny run closure: extracts solver settings from Shiny inputs
    shiny_run = function(model, init_list, params_list, input) {
        run_model(model,
            init = init_list,
            params = params_list,
            times = list(start = input$ode_t0, duration = input$ode_tt, step = input$ode_dt),
            options = list(method = input$ode_method))
    }

    structure(list(
        type = "ODE",
        init = init, params = params, equations = equations, times = times,
        options = modifyList(list(method = "lsoda"), options),
        shiny = list(ui = shiny_ui, run = shiny_run)
    ), class = c("ode_model", "modeller"))
}

#' @export
run_model.ode_model = function(model, init = NULL, params = NULL,
    times = NULL, options = NULL, ...)
{
    init = modifyList(model$init, init %||% list())
    params = modifyList(model$params, params %||% list())
    times = resolve_times(model$times, times %||% list(), has_step = TRUE)
    options = modifyList(model$options, options %||% list())

    tval = seq(times$start, times$stop, times$step)

    sol = deSolve::ode(y = unlist(init), times = tval,
        func = model$equations, parms = params, method = options$method)
    data = as.data.frame(sol)
    names(data)[1] = "t"

    attr(data, "dt") = times$step
    attr(data, "geom") = "line"
    class(data) = c("model_result", class(data))

    standard_checks(data)

    data
}
