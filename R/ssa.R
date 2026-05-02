#' Create a stochastic SSA model
#'
#' Defines a compartmental model solved with Gillespie's stochastic simulation
#' algorithm via [adaptivetau]. Transitions are defined using the `tx()`
#' syntax: `tx(A -> B) = rate` specifies that compartment `A` loses one
#' individual and compartment `B` gains one, at the given rate. For example,
#' `tx(S -> I) = beta * S * I / N`.
#'
#' The equations function is written as a plain R function with no arguments.
#' Inside the function body, compartment names (from `init`) and parameter
#' names (from `params`) can be used as ordinary variables — the package
#' makes them available automatically. The variable `t` (the current time)
#' is also available.
#'
#' @param init Named list of initial compartment values, typically integers
#'   (e.g. `list(S = 999, I = 1, R = 0)`).
#' @param params Named list of parameter values. Must include `time`, which
#'   controls the simulation time span and can be specified as:
#'   - A single number `N` (duration): simulates from 0 to N.
#'   - `c(start, stop)`: simulates from start to stop.
#'   - `c(start, stop, step)`: the step element is stored but ignored by
#'     the SSA solver (stochastic events occur at irregular times).
#'
#'   Inside the equations function, `time` is a three-element numeric vector
#'   `c(start, stop, step)` regardless of how it was originally specified.
#' @param equations A function with no arguments. Use `tx(A -> B) = rate` to
#'   define transitions between compartments. Compartments whose names start
#'   with `total_` are treated as cumulative counters; the corresponding
#'   incidence is computed automatically when plotting.
#' @param options Named list of solver options. For SSA models:
#'   - `method`: `"adaptive-tau"` (default) or `"exact"`.
#'   - `seed`: random number seed, or `NULL` (default) for
#'     non-deterministic results.
#'   - `epsilon`: tau-leaping error control parameter passed to
#'     [adaptivetau::ssa.adaptivetau()] as `tl.params$epsilon`. Smaller
#'     values give more accurate but slower simulations. Default `0.05`.
#'     Ignored when `method = "exact"`.
#'
#'   Default: `list(seed = NULL, method = "adaptive-tau", epsilon = 0.05)`.
#'   Can be overridden in [run_model()].
#'
#' @return An object of class `c("ssa_model", "modeller")`. Use
#'   [run_model()] to simulate, [plot()][plot.model_result] to visualise,
#'   and [show_model()] to launch the interactive app.
#'
#' @examples
#' m = ssa_model(
#'     init = list(S = 999, I = 1, R = 0),
#'     params = list(beta = 0.3, gamma = 0.1, time = 100),
#'     equations = function() {
#'         N = S + I + R
#'         tx(S -> I) = beta * I / N * S
#'         tx(I -> R) = gamma * I
#'     }
#' )
#' m
#' result = run_model(m)
#' plot(result)
#'
#' @export
ssa_model = function(init, params, equations, options = list())
{
    params = validate_inputs(init, params, equations)

    # Build transitions and edit body of equations function
    eq = body(equations)
    txs = elixir::expr_match(eq, { tx(..A -> ..B) = ..C } ? { tx(..A -> ..B) <- ..C })
    transitions = list()

    # Parse a transition side (e.g. A, A + B, or 0) into compartment names.
    parse_compartments = function(expr) {
        if (is.call(expr) && identical(expr[[1]], as.name("+"))) {
            c(parse_compartments(expr[[2]]), parse_compartments(expr[[3]]))
        } else if (is.name(expr)) {
            as.character(expr)
        } else if (is.numeric(expr) && expr == 0) {
            character(0)
        } else {
            stop("malformed transition: expected compartment name or 0, got `",
                deparse(expr), "`", call. = FALSE)
        }
    }

    # Schema for vector compartments. Scalar -> length 1. Each tx() group
    # generates K transitions where K is the length of the vector
    # compartments referenced; scalar compartments broadcast to every
    # transition in the group.
    schema = vapply(init, length, integer(1))
    flat_names = build_flat_names(schema)
    n_flat = length(flat_names)
    flat_offsets = setNames(c(0L, head(cumsum(schema), -1L)), names(schema))

    tvec_offset = 0L
    for (i in seq_along(txs)) {
        tx = txs[[i]]
        from = parse_compartments(tx$A)
        to = parse_compartments(tx$B)

        for (name in c(from, to)) {
            if (!name %in% names(init))
                stop("unknown compartment in transition: `", name, "`", call. = FALSE)
        }

        # Determine K: the length of any vector compartment in this tx.
        # All vector compartments referenced must agree on length.
        cmpts = c(from, to)
        cmpt_lens = unique(schema[cmpts])
        cmpt_lens = cmpt_lens[cmpt_lens > 1]
        if (length(cmpt_lens) > 1) {
            stop("transition references vector compartments of different ",
                "lengths: ", paste(cmpts, collapse = ", "), call. = FALSE)
        }
        K = if (length(cmpt_lens) == 0) 1L else cmpt_lens[[1]]

        # Build K transitions, each a length-n_flat increment vector.
        for (k in seq_len(K)) {
            increments = structure(rep(0, n_flat), names = flat_names)
            for (name in from) {
                idx = if (schema[[name]] == 1L) 1L else k
                pos = flat_offsets[[name]] + idx
                increments[pos] = increments[pos] - 1
            }
            for (name in to) {
                idx = if (schema[[name]] == 1L) 1L else k
                pos = flat_offsets[[name]] + idx
                increments[pos] = increments[pos] + 1
            }
            transitions[[tvec_offset + k]] = increments
        }

        # Rewrite the rate expression. Single transition -> single slot;
        # vector group -> slice of K slots filled by the (length-K) rate
        # expression (or scalar broadcast).
        if (K == 1L) {
            eq[[tx$loc]] = rlang::expr(.tvec[[!!(tvec_offset + 1L)]] <<- !!tx$C)
        } else {
            slot_seq = (tvec_offset + 1L):(tvec_offset + K)
            eq[[tx$loc]] = rlang::expr(.tvec[!!slot_seq] <<- !!tx$C)
        }

        tvec_offset = tvec_offset + K
    }
    check_helpers(eq, "SSA")

    # Put equations function into correct form
    formals(equations) = alist(.state =, .params =, t =)
    body(equations) = rlang::expr({
        record = modeller:::null_record
        .tvec = numeric(!!length(transitions))
        t = t + .params$time[1]
        state_list = modeller:::unpack_state(.state, !!schema)
        with(c(state_list, .params), !!eq)
        return (.tvec)
    })

    # Build recorder function
    recorder = build_recorder(equations)

    # Shiny UI elements for solver settings
    ssa_methods = c("adaptive-tau", "exact")
    default_options = modifyList(list(seed = NULL, method = "adaptive-tau", epsilon = 0.05), options)
    shiny_ui = list(
        shiny::tags$p(shiny::tags$strong("SSA settings")),
        inshiny::inline("Random seed: ", inshiny::inline_number("ssa_seed",
            value = default_options$seed %||% 1,
            min = 1, step = 1, meaning = "Random number generator seed")),
        inshiny::inline("Method: ", inshiny::inline_select("ssa_method",
            choices = ssa_methods, selected = default_options$method,
            meaning = "SSA method")),
        inshiny::inline("Epsilon: ", inshiny::inline_number("ssa_epsilon",
            value = default_options$epsilon, min = 0.01, max = 0.99, step = 0.01,
            arrows = FALSE,
            meaning = "Tau-leaping error control (smaller = more accurate, slower)")),
        inshiny::inline("Start time: ",
            inshiny::inline_number("ssa_t0", value = params$time[1], min = 0,
                placeholder = params$time[1], arrows = FALSE, meaning = "Start time")),
        inshiny::inline("Duration: ",
            inshiny::inline_number("ssa_tt", value = params$time[2] - params$time[1], min = 0,
                placeholder = params$time[2] - params$time[1], arrows = FALSE, meaning = "Duration"))
    )

    # Shiny run closure: extracts solver settings from Shiny inputs
    shiny_run = function(model, init_list, params_list, input) {
        params_list$time = c(input$ssa_t0, input$ssa_t0 + input$ssa_tt, 1)
        run_model(model,
            init = init_list,
            params = params_list,
            options = list(seed = input$ssa_seed, method = input$ssa_method,
                epsilon = input$ssa_epsilon))
    }

    structure(list(
        type = "SSA",
        init = init, params = params,
        equations = equations, recorder = recorder,
        transitions = transitions,
        options = default_options,
        shiny = list(ui = shiny_ui, run = shiny_run)
    ), class = c("ssa_model", "modeller"))
}

#' @export
run_model.ssa_model = function(model, init = NULL, params = NULL,
    options = NULL, ...)
{
    init = modifyList(model$init, init %||% list())
    params = modifyList(model$params, params %||% list())
    params = validate_inputs(init, params, ".BYPASS")
    options = modifyList(model$options, options %||% list())

    if (!is.null(options$seed)) set.seed(options$seed)
    params_list = params

    if (options$method == "adaptive-tau") {
        sol = adaptivetau::ssa.adaptivetau(pack_state(init),
            model$transitions, model$equations, params_list, params$time[2] - params$time[1],
            tl.params = list(epsilon = options$epsilon))
    } else {
        sol = adaptivetau::ssa.exact(pack_state(init),
            model$transitions, model$equations, params_list, params$time[2] - params$time[1])
    }
    data = as.data.frame(sol)
    names(data)[1] = "t"
    data[[1]] = data[[1]] + params$time[1]

    attr(data, "dt") = 0
    attr(data, "geom") = "step"
    data = compute_incidence(data)
    data = compute_recordings(model, params, data)
    data = remove_totals(data)
    class(data) = c("model_result", class(data))

    standard_checks(data)

    data
}


grid_rows = function(x, times)
{
    rows = findInterval(times, x[, 1], rightmost.closed = FALSE)
    ran = range(x[, 1])
    rows[times < ran[1] | times > ran[2]] = NA
    return (rows)
}

grid_summarize = function(x, times, ...)
{
    # Check type of x and force to list of matrices
    if (is.matrix(x)) {
        x = list(x)
        given_list = FALSE
    } else if (is.list(x) && all(sapply(x, is.matrix))) {
        given_list = TRUE
    } else {
        stop("expected matrix or list of matrices for x")
    }

    # Get functions to apply
    funcs = list(...)
    if (given_list && length(funcs) == 0) {
        stop("supply at least one summary function in `...`")
    } else if (!given_list && length(funcs) != 0) {
        stop("do not supply `...` when x is a matrix")
    }

    if (length(funcs) == 0) {
        funcs = list(function(x) x[1])
    }

    for (i in seq_along(funcs)) {
        if (!is.function(funcs[[i]])) {
            stop("expecting a function for each element of `...`")
        }
    }

    # Get rows corresponding to times for each of x
    rows = lapply(x, grid_rows, times = times)

    # Create full 3d matrix
    # x[i, j, k] is now the ith row (for times[i]), jth column, run k.
    x = mapply(function(i, r) x[[i]][r, , drop = FALSE], seq_along(x), rows, SIMPLIFY = "array")

    # Create output matrix
    dcol = ncol(x) - 1 # data columns, same for each run
    tcol = 1 + dcol * length(funcs) # output columns, first of which is time
    y = matrix(times, ncol = tcol, nrow = length(times))
    colnames(y) = rep("time", tcol)

    # Populate matrix
    i = 1
    for (f in seq_along(funcs)) {
        func = funcs[[f]]
        for (c in seq_len(dcol)) {
            i = i + 1
            y[, i] = apply(x[, c + 1, , drop = FALSE], 1, func)
            colnames(y)[i] = paste0(colnames(x)[c + 1], names(funcs)[f])
        }
    }

    return (y)
}
