#' Create a discrete-time difference equation model
#'
#' Defines a compartmental model solved by iterating a difference equation.
#' At each time step, the equations function calculates and returns the
#' next value of every compartment. For example, if compartment `S` starts
#' at 1000 and 10 individuals are infected this step, the equations should
#' return `S = S - 10` (i.e. 990), not the change (`-10`).
#'
#' The equations function is written as a plain R function with no arguments.
#' Inside the function body, compartment names (from `init`) and parameter
#' names (from `params`) can be used as ordinary variables — the package
#' makes them available automatically. The variable `dt` (the time step
#' size) is also available.
#'
#' @param init Named list of initial compartment values (e.g.
#'   `list(S = 999, I = 1, R = 0)`).
#' @param params Named list of parameter values (e.g.
#'   `list(beta = 0.3, gamma = 0.1)`).
#' @param equations A function with no arguments. Inside the function body,
#'   write equations using compartment and parameter names as if they were
#'   ordinary variables. Return a list giving the next value of each
#'   compartment, in the same order as `init`. Compartments whose names start
#'   with `total_` are treated as cumulative counters; the corresponding
#'   incidence is computed automatically when plotting.
#' @param times Named list controlling the time span. Accepts `start`,
#'   `stop`, `duration`, and `step`. At minimum, provide `duration` or
#'   `stop`. Defaults: `start = 0`, `stop = 100`, `step = 1`.
#' @param options Named list of solver options. No solver options are
#'   available for difference equation models.
#'
#' @return An object of class `c("difference_model", "modeller")`. Use
#'   [run_model()] to simulate, [plot()][plot.model_result] to visualise,
#'   and [show_model()] to launch the interactive app.
#'
#' @examples
#' m = difference_model(
#'     init = list(S = 999, I = 1, R = 0),
#'     params = list(beta = 0.3, gamma = 0.1),
#'     equations = function() {
#'         N = S + I + R
#'         lambda = beta * I / N
#'         list(
#'             S = S - lambda * S,
#'             I = I + lambda * S - gamma * I,
#'             R = R + gamma * I
#'         )
#'     },
#'     times = list(duration = 100)
#' )
#' m
#' result = run_model(m)
#' plot(result)
#'
#' @export
difference_model = function(init, params, equations, times, options = list())
{
    validate_inputs(init, params, equations, times)

    # Put equations function into correct form
    formals(equations) = alist(t =, .state =, .params =)
    eq = body(equations)
    eq = elixir::expr_replace(eq, { new(`.A:name`) <- ..B },
        { .nlist[[as.character(quote(.A))]] <<- ..B })
    eq = elixir::expr_replace(eq, { new(`.A:name`) = ..B },
        { .nlist[[as.character(quote(.A))]] <<- ..B })

    body(equations) = rlang::expr({
        .nlist = init
        with(c(list(dt = .params$.dt), .state, .params), !!eq)
        return (unlist(.nlist))
    })

    # Process times
    times = resolve_times(times, has_step = TRUE)

    # Shiny UI elements for solver settings
    shiny_ui = list(
        shiny::tags$p(shiny::tags$strong("Difference equation settings")),
        inshiny::inline("Start time: ",
            inshiny::inline_number("diff_t0", value = times$start, min = 0,
                placeholder = times$start, arrows = FALSE, meaning = "Start time")),
        inshiny::inline("Duration: ",
            inshiny::inline_number("diff_tt", value = times$duration, min = 0,
                placeholder = times$duration, arrows = FALSE, meaning = "Duration")),
        inshiny::inline("Time step: ",
            inshiny::inline_number("diff_dt", value = times$step, min = 0.01, step = 0.01,
                placeholder = times$step, arrows = FALSE, meaning = "Time step"))
    )

    # Shiny run closure: extracts solver settings from Shiny inputs
    shiny_run = function(model, init_list, params_list, input) {
        run_model(model,
            init = init_list,
            params = params_list,
            times = list(start = input$diff_t0, duration = input$diff_tt, step = input$diff_dt))
    }

    structure(list(
        type = "difference",
        init = init, params = params, equations = equations, times = times,
        options = options,
        shiny = list(ui = shiny_ui, run = shiny_run)
    ), class = c("difference_model", "modeller"))
}

#' @export
run_model.difference_model = function(model, init = NULL, params = NULL,
    times = NULL, options = NULL, ...)
{
    init = modifyList(model$init, init %||% list())
    params = modifyList(model$params, params %||% list())
    times = resolve_times(model$times, times %||% list(), has_step = TRUE)

    tval = seq(times$start, times$stop, times$step)
    state = unlist(init)
    params_list = c(params, list(.dt = times$step))
    n = length(tval)

    # Pre-allocate matrix: rows = time steps, cols = t + compartments
    result = matrix(NA_real_, nrow = n, ncol = 1L + length(state))
    result[1, ] = c(tval[1], state)

    for (i in 2:n) {
        state = model$equations(tval[i - 1], state, params_list)
        result[i, ] = c(tval[i], state)
    }

    data = as.data.frame(result)
    names(data) = c("t", names(init))

    attr(data, "dt") = times$step
    attr(data, "geom") = "step"
    class(data) = c("model_result", class(data))

    data = compute_incidence(data)
    standard_checks(data)

    data
}
