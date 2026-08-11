#' Create a deterministic ODE model
#'
#' Defines a compartmental model solved with ordinary differential equations
#' via [deSolve::ode()]. At each time point, the equations function
#' calculates the rate of change (derivative) of every compartment using the
#' `d()` syntax. For example, if individuals leave compartment `S` at rate
#' `alpha * S`, write `d(S) = -alpha * S`.
#'
#' The equations function is written as a plain R function with no arguments.
#' Inside the function body, compartment names (from `init`) and parameter
#' names (from `params`) can be used as ordinary variables — the package
#' makes them available automatically. The variable `t` (the current time)
#' is also available.
#'
#' @param init Named list of initial compartment values (e.g.
#'   `list(S = 999, I = 1, R = 0)`).
#' @param params Named list of parameter values. Must include `time`, which
#'   controls the simulation time span and can be specified as:
#'   - A single number `N` (duration): simulates from 0 to N with step 1.
#'   - `c(start, stop)`: simulates from start to stop with step 1.
#'   - `c(start, stop, step)`: simulates from start to stop with the given
#'     step size.
#'
#'   Inside the equations function, `time` is a three-element numeric vector
#'   `c(start, stop, step)` regardless of how it was originally specified.
#' @param equations A function with no arguments. Use `d(X) = ...` to set
#'   the derivative (rate of change) of each compartment `X`. Compartments
#'   whose names start with `total_` are treated as cumulative counters; the
#'   corresponding incidence is computed automatically when plotting.
#' @param options Named list of solver options. For ODE models the only
#'   option is `method`, the integration method passed to [deSolve::ode()].
#'   Default: `list(method = "lsode")`. Can be overridden in [run_model()].
#'
#' @return An object of class `c("ode_model", "modeller")`. Use
#'   [run_model()] to simulate, [plot()][plot.model_result] to visualise,
#'   and [show_model()] to launch the interactive app.
#'
#' @examples
#' m = ode_model(
#'     init = list(S = 999, I = 1, R = 0),
#'     params = list(beta = 0.3, gamma = 0.1, time = 100),
#'     equations = function() {
#'         N = S + I + R
#'         d(S) = -beta * I / N * S
#'         d(I) =  beta * I / N * S - gamma * I
#'         d(R) =  gamma * I
#'     }
#' )
#' m
#' result = run_model(m)
#' plot(result)
#'
#' @export
ode_model = function(init, params, equations, options = list())
{
    params = validate_inputs(init, params, equations)

    # Put equations function into correct form
    formals(equations) = alist(t =, .state =, .params =)
    eq = body(equations)

    # Catch d() calls referring to unknown compartments
    d_calls = elixir::expr_match(eq,
        { d(`.A:name`) = ..B } ? { d(`.A:name`) <- ..B })
    d_names = unique(vapply(d_calls, function(m) deparse(m$A), character(1)))
    unknown = setdiff(d_names, names(init))
    if (length(unknown) > 0) {
        stop("d() refers to unknown compartment(s): ",
            paste(unknown, collapse = ", "), call. = FALSE)
    }

    eq = elixir::expr_replace(eq, { d(`.A:name`) <- ..B },
        { .dlist[[as.character(quote(.A))]] <<- ..B })
    eq = elixir::expr_replace(eq, { d(`.A:name`) = ..B },
        { .dlist[[as.character(quote(.A))]] <<- ..B })
    check_helpers(eq, "ODE")

    # Name unnamed record() arguments.
    eq = name_record_args(eq)

    # Schema for packing/unpacking vector compartments. Scalar compartments
    # (length 1) are unaffected; vector ones flatten into Name.1..Name.K.
    schema = vapply(init, length, integer(1))

    # .dlist defaults to zero-vectors of the right length so a compartment
    # without d() stays constant (rate of change = 0).
    unpack_expr = build_unpack_expr(schema)

    flat_names = build_flat_names(schema)
    body(equations) = rlang::expr({
        record = function(...) NULL
        .dlist = !!lapply(init, function(x) numeric(length(x)))
        state_list = !!unpack_expr
        with(c(state_list, .params), !!eq)
        .flat = unlist(.dlist, use.names = FALSE)
        names(.flat) = !!flat_names
        return (list(.flat))
    })

    # Build recorder function
    recorder = build_recorder(equations)

    # Shiny UI elements for solver settings
    ode_methods = c("lsoda", "lsode", "lsodes", "lsodar", "vode", "daspk",
           "euler", "rk4", "ode23", "ode45", "radau",
           "bdf", "bdf_d", "adams", "impAdams", "impAdams_d")
    default_options = modifyList(list(method = "lsode"), options)
    shiny_ui = list(
        shiny::tags$p(shiny::tags$strong("ODE settings")),
        inshiny::inline("Method: ", inshiny::inline_select("ode_method",
            choices = ode_methods, selected = default_options$method,
            meaning = "ODE method")),
        inshiny::inline("Start time: ",
            inshiny::inline_number("ode_t0", value = params$time[1], min = 0,
                placeholder = params$time[1], arrows = FALSE, meaning = "Start time")),
        inshiny::inline("Duration: ",
            inshiny::inline_number("ode_tt", value = params$time[2] - params$time[1], min = 0,
                placeholder = params$time[2] - params$time[1], arrows = FALSE, meaning = "Duration")),
        inshiny::inline("Time step: ",
            inshiny::inline_number("ode_dt", value = params$time[3], min = 0.01, step = 0.01,
                placeholder = params$time[3], arrows = FALSE, meaning = "Time step"))
    )

    # Shiny run closure: extracts solver settings from Shiny inputs
    shiny_run = function(model, init_list, params_list, input) {
        params_list$time = c(input$ode_t0, input$ode_t0 + input$ode_tt, input$ode_dt)
        run_model(model,
            init = init_list,
            params = params_list,
            options = list(method = input$ode_method))
    }

    # Default values for solver-settings inputs (used by the Reset button)
    shiny_defaults = list(
        ode_method = default_options$method,
        ode_t0 = params$time[1],
        ode_tt = params$time[2] - params$time[1],
        ode_dt = params$time[3]
    )

    structure(list(
        type = "ODE",
        init = init, params = params,
        equations = equations, recorder = recorder,
        options = modifyList(list(method = "lsode"), options),
        shiny = list(ui = shiny_ui, run = shiny_run, defaults = shiny_defaults)
    ), class = c("ode_model", "modeller"))
}

#' @export
run_model.ode_model = function(model, init = NULL, params = NULL,
    options = NULL, ...)
{
    init = modifyList(model$init, init %||% list())
    params = modifyList(model$params, params %||% list())
    params = validate_inputs(init, params, ".BYPASS")
    options = modifyList(model$options, options %||% list())

    tval = seq(params$time[1], params$time[2], params$time[3])
    if (tval[length(tval)] < params$time[2]) tval = c(tval, params$time[2])

    sol = deSolve::ode(y = pack_state(init), times = tval,
        func = model$equations, parms = params, method = options$method)
    data = as.data.frame(sol)
    names(data)[1] = "t"

    attr(data, "dt") = params$time[3]
    attr(data, "geom") = "line"
    data = compute_incidence(data)
    data = compute_recordings(model, params, data)
    data = remove_totals(data)
    class(data) = c("model_result", class(data))

    standard_checks(data)

    data
}
