#' @export
ssa_model = function(init, params, transitions, model, times, ...)
{
    # Put model function into correct form
    # TODO check function's signature
    formals(model) = alist(.state =, .params =, t =)
    body(model) = rlang::expr({
        t = t + .params$.tstart
        .x = with(as.list(c(.state, .params)), !!body(model))
        return (unlist(.x))
    })

    # Process times
    times$start = times$start %||% 0
    times$stop = times$stop %||% 100
    times$duration = times$duration %||% (times$stop - times$start)

    # Extra elements for model configuration
    ssa_methods = c("adaptive-tau", "exact")
    extra_elements = list(
        shiny::h5("SSA settings"),
        inshiny::inline("Random seed: ", inshiny::inline_number("ssa_seed", 1,
            min = 1, step = 1, meaning = "Random number generator seed")),
        inshiny::inline("Method: ", inshiny::inline_select("ssa_method",
            choices = ssa_methods, meaning = "SSA method")),
        inshiny::inline("Start time: ",
            inshiny::inline_number("ssa_t0", value = times$start, min = 0,
                placeholder = times$start, arrows = FALSE, meaning = "Start time")),
        inshiny::inline("Duration: ",
            inshiny::inline_number("ssa_tt", value = times$duration, min = 0,
                placeholder = times$duration, arrows = FALSE, meaning = "Duration"))
    )

    # Solver: returns data frame with t + compartment columns
    run_model = function(init_list, params_list, times, input) {
        set.seed(input$ssa_seed)
        params_list = c(params_list, list(.tstart = input$ssa_t0))
        if (input$ssa_method == "adaptive-tau") {
            sol = adaptivetau::ssa.adaptivetau(unlist(init_list), transitions,
                model, params_list, input$ssa_tt)
        } else {
            sol = adaptivetau::ssa.exact(unlist(init_list), transitions,
                model, params_list, input$ssa_tt)
        }
        data = as.data.frame(sol)
        names(data)[1] = "t"
        data[[1]] = data[[1]] + input$ssa_t0
        attr(data, "dt") = 0
        data
    }

    model_app(run_model, init, params, times, extra_elements, ...)
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
