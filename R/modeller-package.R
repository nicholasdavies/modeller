#' Interactive dynamical systems modelling
#'
#' Build, run, and interactively explore compartmental dynamical systems
#' models using deterministic (ODE), stochastic (SSA), and discrete-time
#' (difference equation) solvers.
#'
#' The typical workflow is:
#'
#' 1. Define a model with [ode_model()], [ssa_model()], or [difference_model()]
#' 2. Run it with [run_model()]
#' 3. Visualise results with [plot()][plot.model_result]
#' 4. Explore interactively with [show_model()]
#'
#' @section Model specification:
#' Equations are written as specially-formatted R functions that reference
#' compartment and parameter names directly. The package transforms these into
#' the form needed by the underlying solvers (`deSolve` for ODEs, `adaptivetau`
#' for SSA, simple iteration for difference equations).
#'
#' @keywords internal
#' @importFrom utils head modifyList
#' @importFrom stats optim plogis qlogis
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL

# Symbols used in non-standard evaluation contexts (elixir pattern matching
# in expr_match/expr_replace and ggplot2 aesthetic mappings) that R CMD
# check otherwise flags as undefined globals.
utils::globalVariables(c(
    "?", "d", "d<-", "new", "new<-", "tx",
    "..A", "..B", "...A", ".params",
    ".data", "name", "value"
))

#' Run a model
#'
#' Runs a model with default or overridden parameters, initial conditions,
#' and time settings. Returns a data frame of class `model_result`.
#'
#' Any argument except `model` can be omitted to use the defaults stored in
#' the model object. Partial overrides are supported: e.g.
#' `params = list(R0 = 5)` overrides only `R0`, keeping all other parameters
#' at their defaults.
#'
#' @param model A model object created by [ode_model()], [ssa_model()], or
#'   [difference_model()].
#' @param init Named list of initial conditions to override defaults.
#' @param params Named list of parameters to override defaults.
#' @param options Named list of solver options to override defaults.
#'
#'   For ODE models (default: `list(method = "lsode")`):
#'   - `method`: integration method passed to [deSolve::ode()]. Common
#'     choices include `"lsode"` (default, adaptive Adams/BDF), `"lsoda"`
#'     (adaptive with automatic stiffness switching), `"rk4"` (fixed-step
#'     Runge-Kutta), and `"euler"`.
#'
#'   For SSA models (default: `list(seed = NULL, method = "adaptive-tau",
#'   epsilon = 0.05)`):
#'   - `seed`: random number seed for reproducibility, or `NULL` (default)
#'     for non-deterministic results.
#'   - `method`: `"adaptive-tau"` (default, faster approximate) or `"exact"`
#'     (Gillespie's direct method).
#'   - `epsilon`: tau-leaping error control. Smaller is more accurate but
#'     slower. Default `0.05`. Ignored when `method = "exact"`.
#'
#'   For difference equation models: no solver options are available.
#' @param ... Ignored.
#'
#' @return A `model_result` data frame with a `t` column and one column per
#'   compartment. Has attributes `dt` (time step; 0 for SSA) and `geom`
#'   (`"line"` for ODE, `"step"` for SSA and difference equation models).
#'
#' @examples
#' m = ode_model(
#'     init = list(S = 999, I = 1),
#'     params = list(beta = 0.3, gamma = 0.1, time = 50),
#'     equations = function() {
#'         d(S) = -beta * S * I
#'         d(I) =  beta * S * I - gamma * I
#'     }
#' )
#' result = run_model(m)
#' result = run_model(m, params = list(beta = 0.5))
#' result = run_model(m, options = list(method = "rk4"))
#'
#' @export
run_model = function(model, init = NULL, params = NULL,
    options = NULL, ...) UseMethod("run_model")

#' @export
print.modeller = function(x, ...)
{
    type = if (!is.null(x$type)) x$type else "Unknown"

    param_names = setdiff(names(x$params), "time")

    fmt_compartment = function(nm, val) {
        if (length(val) == 1) sprintf("%s(0) = %s", nm, format(val))
        else if (all(val == val[1])) sprintf("%s[%d](0) = %s", nm, length(val), format(val[1]))
        else sprintf("%s[%d](0) = c(%s, ...)", nm, length(val), format(val[1]))
    }
    fmt_param = function(nm, val) {
        if (length(val) == 1) sprintf("%s = %s", nm, format(val))
        else sprintf("%s[%d] = c(%s, ...)", nm, length(val), format(val[1]))
    }

    cat(type, "model\n")
    cat("  Compartments: ",
        paste(mapply(fmt_compartment, names(x$init), x$init), collapse = ", "),
        "\n", sep = "")
    cat("  Parameters:   ",
        paste(mapply(fmt_param, param_names, x$params[param_names]), collapse = ", "),
        "\n", sep = "")
    cat("  Time:         ", x$params$time[1], " to ", x$params$time[2], ", step ", x$params$time[3], "\n", sep = "")
    invisible(x)
}

# Check that the equations body doesn't reference helpers from other model
# types. Called from each model constructor after that model's own helper
# rewrites have been applied; any remaining d()/new()/tx() calls indicate
# either the wrong helper for this model type or a malformed call.
check_helpers = function(eq, model_type)
{
    if (model_type != "ODE" && elixir::expr_detect(eq, { d(..A) })) {
        stop("`d()` is the syntax for ODE models. Use ",
            if (model_type == "SSA") "`tx()`" else "`new()`",
            " for ", model_type, " models.", call. = FALSE)
    }
    if (model_type != "difference" && elixir::expr_detect(eq, { new(..A) })) {
        stop("`new()` is the syntax for difference equation models. Use ",
            if (model_type == "SSA") "`tx()`" else "`d()`",
            " for ", model_type, " models.", call. = FALSE)
    }
    if (elixir::expr_detect(eq, { tx(..A) })) {
        if (model_type != "SSA") {
            stop("`tx()` is the syntax for SSA models. Use ",
                if (model_type == "ODE") "`d()`" else "`new()`",
                " for ", model_type, " models.", call. = FALSE)
        } else {
            stop("malformed `tx()` call. Expected `tx(A -> B) = rate`.",
                call. = FALSE)
        }
    }
}

# Validate inputs common to all model constructors
validate_inputs = function(init, params, equations)
{
    # Basic structure of arguments
    if (!is.list(init) || is.null(names(init)) || any(names(init) == ""))
        stop("`init` must be a named list", call. = FALSE)
    if (!all(vapply(init, is.numeric, logical(1))))
        stop("all elements of `init` must be numeric", call. = FALSE)
    if (!is.list(params) || is.null(names(params)) || any(names(params) == ""))
        stop("`params` must be a named list", call. = FALSE)
    if (!all(vapply(params, is.numeric, logical(1))))
        stop("all elements of `params` must be numeric", call. = FALSE)
    if (!identical(equations, ".BYPASS")) {
        if (!is.function(equations))
            stop("`equations` must be a function", call. = FALSE)
        if (length(formals(equations)) > 0)
            stop("`equations` must be a function with no arguments", call. = FALSE)
    }

    # Reserved names: `t` is the simulation time, `dt` is the difference
    # equation step, names starting with `.` are internal.
    reserved = c("t", "dt")
    is_internal = function(n) startsWith(n, ".")

    bad_init = names(init) %in% c(reserved, "time") | vapply(names(init), is_internal, logical(1))
    if (any(bad_init)) {
        stop("`init` names cannot be `t`, `dt`, `time`, or start with `.`: ",
            paste(names(init)[bad_init], collapse = ", "), call. = FALSE)
    }
    bad_params = names(params) %in% reserved | vapply(names(params), is_internal, logical(1))
    if (any(bad_params)) {
        stop("`params` names cannot be `t`, `dt`, or start with `.`: ",
            paste(names(params)[bad_params], collapse = ", "), call. = FALSE)
    }
    overlap = intersect(names(init), names(params))
    if (length(overlap) > 0) {
        stop("`init` and `params` cannot share names: ",
            paste(overlap, collapse = ", "), call. = FALSE)
    }

    # Times
    if (length(params$time) == 0) {
        warning("no `time` parameter supplied; setting duration to 100.")
        params$time = c(0, 100, 1)
    } else if (length(params$time) == 1 && params$time > 0) {
        params$time = c(0, params$time, min(params$time, 1))
    } else if (length(params$time) == 2 && params$time[2] > params$time[1]) {
        params$time = c(params$time, min(params$time[2] - params$time[1], 1))
    } else if (length(params$time) == 3 && params$time[2] > params$time[1] && params$time[3] > 0) {
        params$time = c(params$time[1], params$time[2], min(params$time[2] - params$time[1], params$time[3]))
    } else {
        stop("`params$time` must either be: duration; c(start_time, stop_time); or c(start_time, stop_time, time_step)", call. = FALSE)
    }

    return (params)
}

# Build the flat state-vector column names from a schema (named integer
# vector of compartment lengths). Scalars get the bare name; vectors get
# "name.1", "name.2", ... (R 1-indexed for consistency with R subscripting).
build_flat_names = function(schema)
{
    parts = character(0)
    for (i in seq_along(schema)) {
        nm = names(schema)[i]
        n = schema[[i]]
        if (n == 1L) {
            parts = c(parts, nm)
        } else {
            parts = c(parts, paste0(nm, ".", seq_len(n)))
        }
    }
    parts
}

# Pack a named list of vectors into a single named numeric vector, using
# the build_flat_names convention. The complement of unpack_state.
pack_state = function(state_list)
{
    schema = vapply(state_list, length, integer(1))
    flat = unlist(state_list, use.names = FALSE)
    names(flat) = build_flat_names(schema)
    flat
}

# Unpack a flat numeric vector into a named list of vectors using a schema.
# The complement of pack_state.
# Build a literal list expression that slices a flat state vector into a
# named list of compartments, with the slice indices baked in at model
# construction time. Used in equations bodies to avoid the per-call cost of
# unpack_state's schema loop.
build_unpack_expr = function(schema)
{
    ends = cumsum(schema)
    starts = c(1L, head(ends, -1L) + 1L)
    parts = sprintf("%s = .state[%d:%d]", names(schema), starts, ends)
    parse(text = sprintf("list(%s)", paste(parts, collapse = ", ")))[[1]]
}

unpack_state = function(flat, schema)
{
    out = vector("list", length(schema))
    names(out) = names(schema)
    idx = 1L
    for (i in seq_along(schema)) {
        n = schema[[i]]
        v = flat[idx:(idx + n - 1L)]
        names(v) = NULL
        out[[i]] = v
        idx = idx + n
    }
    out
}

# Prefix for cumulative counter compartments that get converted to incidence.
cumulative_prefix = "total_"

# Convert cumulative total_ columns to incidence rates. Removes total_
# columns and adds corresponding incidence columns (without the prefix).
# For SSA results (irregular time points), computes incidence on a regular
# grid and merges back, using last observation carried forward for the
# original compartment columns at grid time points.
compute_incidence = function(data)
{
    prefix = paste0("^", cumulative_prefix)
    total_cols = grep(prefix, names(data), value = TRUE)
    if (length(total_cols) == 0) return(data)

    dt = attr(data, "dt")
    orig_class = class(data)
    orig_geom = attr(data, "geom")
    grid_step = if (dt <= 1) 1 else dt

    if (dt == 0) {
        # SSA: irregular time points. Compute incidence on a regular grid,
        # then merge with original data and fill NAs with LOCF.
        mat = as.matrix(data[c("t", total_cols)])
        grid = seq(min(data$t), max(data$t), by = grid_step)
        summ = grid_summarize(mat, grid)

        inc_data = data.frame(t = grid)
        for (col in total_cols) {
            new_name = sub(prefix, "", col)
            vals = summ[, col] / grid_step
            inc_data[[new_name]] = c(vals[1], diff(vals))
        }

        data = merge(data, inc_data, by = "t", all = TRUE)
        data = data[order(data$t), ]

        # LOCF for compartment columns at grid-only time points
        for (col in names(data)) {
            if (col != "t" && anyNA(data[[col]])) {
                nna = !is.na(data[[col]])
                data[[col]] = data[[col]][nna][cumsum(nna)]
            }
        }
        rownames(data) = NULL
    } else {
        # Regular grid: diff and divide by dt to get rate per unit time
        for (col in total_cols) {
            new_name = sub(prefix, "", col)
            vals = data[[col]] / dt
            data[[new_name]] = c(vals[1], diff(vals))
        }
    }

    attr(data, "dt") = dt
    attr(data, "geom") = orig_geom
    class(data) = orig_class
    data
}

# Remove totals columns from data
remove_totals = function(data)
{
    prefix = paste0("^", cumulative_prefix)
    total_cols = grep(prefix, names(data), value = TRUE)
    data[total_cols] = NULL
    return (data)
}

# Compute "record" statements for each row of the data.
compute_recordings = function(model, params, data)
{
    if (is.null(model$recorder)) return(data)

    schema = vapply(model$init, length, integer(1))
    flat_state_names = build_flat_names(schema)
    # Other columns (e.g. incidence) are passed via .params so they appear
    # alongside compartments inside the recorder body's `with(...)` env.
    other_cols = setdiff(names(data), c("t", flat_state_names))

    state_mat = as.matrix(data[, flat_state_names, drop = FALSE])
    other_data = lapply(other_cols, function(col) data[[col]])
    names(other_data) = other_cols
    t_vec = data[[1]]
    n = length(t_vec)

    recordings = vector(mode = "list", n)
    if (model$type == "ODE" || model$type == "difference") {
        for (r in seq_len(n)) {
            params_r = params
            for (col in other_cols) params_r[[col]] = other_data[[col]][r]
            recordings[[r]] = model$recorder(t_vec[r], state_mat[r, ], params_r)
        }
    } else if (model$type == "SSA") {
        # SSA equations body converts t from relative to absolute by adding
        # params$time[1]. Recorder shares that body, so we pass relative t
        # here (data t is already absolute) to avoid double-counting.
        start_time = params$time[1]
        for (r in seq_len(n)) {
            params_r = params
            for (col in other_cols) params_r[[col]] = other_data[[col]][r]
            recordings[[r]] = model$recorder(state_mat[r, ], params_r, t_vec[r] - start_time)
        }
    } else {
        stop("Unknown model type.")
    }
    recordings = as.data.frame(data.table::rbindlist(recordings, fill = TRUE))
    result = cbind(data, recordings)
    attr(result, "dt") = attr(data, "dt")
    attr(result, "geom") = attr(data, "geom")
    result
}

# Build recorder function from equations function.
# The first statement of body(equations) is the no-op `record` assignment,
# which is swapped for one that captures values into a local `.rlist`. The
# return statement is then rewritten to return that list.
build_recorder = function(equations)
{
    recorder = NULL
    if (elixir::expr_detect(body(equations), { record(...A) })) {
        recorder = equations
        rec = body(recorder)
        rec[[2]] = rlang::expr({
            .rlist = list()
            record = function(...) .rlist <<- c(.rlist, list(...))
        })
        rec = elixir::expr_replace(rec, { return(..B) }, { return (.rlist) })
        body(recorder) = rec
    }
    return (recorder)
}

# Helper to detect and issue a warning about some problem with the model
# results.
detect_problem = function(result, predicate, description, fix)
{
    # apply() with margin = 2 was the bottleneck for large outputs (it
    # builds a transposed copy via aperm). Coerce once to a matrix and let
    # the predicate work on the whole numeric vector.
    mat = as.matrix(result[, -1, drop = FALSE])
    flagged = predicate(mat)
    if (any(flagged, na.rm = TRUE)) {
        first_row = which(flagged, arr.ind = TRUE)[1, 1]
        warning(description, " detected at time t = ",
            result[first_row, 1], ". ", fix, call. = FALSE)
    }
}

standard_checks = function(result)
{
    detect_problem(result, \(x) x < 0, "Negative compartment",
        "This may be because parameters or initial conditions are negative, step size is too large, or the model logic is incorrect.")
    detect_problem(result, \(x) !is.finite(x), "Non-finite number",
        paste("This may be because the model is generating NAs or conducting undefined mathematical operations (such as log(-1) or division by zero).",
            "Some ODE-solving methods may generate undefined numbers."))
}
