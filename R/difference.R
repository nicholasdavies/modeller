#' Create a discrete-time difference equation model
#'
#' Defines a compartmental model solved by iterating a difference equation.
#' At each time step, the equations function calculates the next value of
#' every compartment using the `new()` syntax. For example, if individuals
#' leave compartment `S` at rate `alpha * S`, write
#' `new(S) = S - alpha * S * dt`.
#'
#' The equations function is written as a plain R function with no arguments.
#' Inside the function body, compartment names (from `init`) and parameter
#' names (from `params`) can be used as ordinary variables — the package
#' makes them available automatically. The variables `t` (the current time)
#' and `dt` (the time step size) are also available.
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
#' @param equations A function with no arguments. Use `new(X) = ...` to set
#'   the next value of each compartment `X`. Compartments whose names start
#'   with `total_` are treated as cumulative counters; the corresponding
#'   incidence is computed automatically when plotting.
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
#'     params = list(beta = 0.3, gamma = 0.1, time = 100),
#'     equations = function() {
#'         N = S + I + R
#'         lambda = beta * I / N
#'         new(S) = S - lambda * S
#'         new(I) = I + lambda * S - gamma * I
#'         new(R) = R + gamma * I
#'     }
#' )
#' m
#' result = run_model(m)
#' plot(result)
#'
#' @export
difference_model = function(init, params, equations, options = list())
{
    # Validate inputs
    params = validate_inputs(init, params, equations)

    # Put equations function into correct form
    formals(equations) = alist(t =, .state =, .params =)
    eq = body(equations)

    # Catch new() calls referring to unknown compartments
    new_calls = elixir::expr_match(eq,
        { new(`.A:name`) = ..B } ? { new(`.A:name`) <- ..B })
    new_names = unique(vapply(new_calls, function(m) deparse(m$A), character(1)))
    unknown = setdiff(new_names, names(init))
    if (length(unknown) > 0) {
        stop("new() refers to unknown compartment(s): ",
            paste(unknown, collapse = ", "), call. = FALSE)
    }

    eq = elixir::expr_replace(eq, { new(`.A:name`) <- ..B },
        { .nlist[[as.character(quote(.A))]] <<- ..B })
    eq = elixir::expr_replace(eq, { new(`.A:name`) = ..B },
        { .nlist[[as.character(quote(.A))]] <<- ..B })
    check_helpers(eq, "difference")

    # Name unnamed record() arguments.
    eq = name_record_args(eq)

    schema = vapply(init, length, integer(1))

    # Pre-compute the unpack as a literal list expression with the state
    # slice indices baked in. Avoids a function call and a per-call schema
    # loop on every iteration.
    unpack_expr = build_unpack_expr(schema)

    # To avoid a second pass for `record()` calls, capture both in first run
    # and keep model$recorder set to NULL.
    body(equations) = rlang::expr({
        .rlist = list()
        record = function(...) .rlist <<- c(.rlist, list(...))
        state_list = !!unpack_expr
        .nlist = state_list
        with(c(list(dt = .params$.dt), state_list, .params), !!eq)
        # use.names = FALSE: set names at end of run_model.difference_model
        list(state = unlist(.nlist, use.names = FALSE), records = .rlist)
    })

    # No separate recorder for difference: records are produced inline.
    recorder = NULL

    # Shiny UI elements for solver settings
    shiny_ui = list(
        shiny::tags$p(shiny::tags$strong("Difference equation settings")),
        inshiny::inline("Start time: ",
            inshiny::inline_number("diff_t0", value = params$time[1], min = 0,
                placeholder = params$time[1], arrows = FALSE, meaning = "Start time")),
        inshiny::inline("Duration: ",
            inshiny::inline_number("diff_tt", value = params$time[2] - params$time[1], min = 0,
                placeholder = params$time[2] - params$time[1], arrows = FALSE, meaning = "Duration")),
        inshiny::inline("Time step: ",
            inshiny::inline_number("diff_dt", value = params$time[3], min = 0.01, step = 0.01,
                placeholder = params$time[3], arrows = FALSE, meaning = "Time step"))
    )

    # Shiny run closure: extracts solver settings from Shiny inputs
    shiny_run = function(model, init_list, params_list, input) {
        params_list$time = c(input$diff_t0, input$diff_t0 + input$diff_tt, input$diff_dt)
        run_model(model,
            init = init_list,
            params = params_list)
    }

    # Default values for solver-settings inputs (used by the Reset button)
    shiny_defaults = list(
        diff_t0 = params$time[1],
        diff_tt = params$time[2] - params$time[1],
        diff_dt = params$time[3]
    )

    structure(list(
        type = "difference",
        init = init, params = params,
        equations = equations, recorder = recorder,
        options = options,
        shiny = list(ui = shiny_ui, run = shiny_run, defaults = shiny_defaults)
    ), class = c("difference_model", "modeller"))
}

#' @export
run_model.difference_model = function(model, init = NULL, params = NULL,
    options = NULL, ...)
{
    init = modifyList(model$init, init %||% list())
    params = modifyList(model$params, params %||% list())
    params = validate_inputs(init, params, ".BYPASS")
    params$.dt = params$time[3]

    tval = seq(params$time[1], params$time[2], params$time[3])
    # Names not needed for the simulation loop; output column names are set
    # at the end. Slicing an unnamed vector is faster than a named one.
    state = unlist(init, use.names = FALSE)
    n = length(tval)

    # Pre-allocate matrix and per-row records list. The equations function
    # returns both next-state and any record() values produced this step.
    # records[[i]] corresponds to the state AT row i (the input to the step,
    # not the output), to match the post-hoc semantics of compute_recordings.
    result = matrix(NA_real_, nrow = n, ncol = 1L + length(state))
    result[, 1] = tval
    cols_state = seq_len(length(state)) + 1L  # state cols 2..end
    result[1, cols_state] = state
    all_records = vector("list", n)

    if (n >= 2) {
        for (i in seq_len(n - 1L)) {
            res = model$equations(tval[i], state, params)
            all_records[[i]] = res$records
            state = res$state
            result[i + 1L, cols_state] = state
        }
    }
    # Final call to capture records for the last output row's state. The
    # returned next-state is discarded (beyond the simulation horizon).
    res = model$equations(tval[n], state, params)
    all_records[[n]] = res$records

    data = as.data.frame(result)
    names(data) = c("t", build_flat_names(vapply(model$init, length, integer(1))))

    attr(data, "dt") = params$time[3]
    attr(data, "geom") = "step"
    data = compute_incidence(data)
    if (any(lengths(all_records) > 0)) {
        records_df = as.data.frame(data.table::rbindlist(all_records, fill = TRUE))
        old_attrs = attributes(data)[c("dt", "geom")]
        data = cbind(data, records_df)
        attr(data, "dt") = old_attrs$dt
        attr(data, "geom") = old_attrs$geom
    }
    data = remove_totals(data)
    class(data) = c("model_result", class(data))

    standard_checks(data)

    data
}
