#' Interactive Infectious Disease Modelling
#'
#' Build, run, and interactively explore compartmental infectious disease
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
#' Equations are written as plain R functions that reference compartment and
#' parameter names directly. The package transforms these into the form needed
#' by the underlying solvers ([deSolve] for ODEs, [adaptivetau] for SSA,
#' simple iteration for difference equations).
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL

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
#' @param times Named list with `start`, `stop`, `duration`, and/or `step`
#'   to override defaults. If `duration` is provided, `stop` is recomputed;
#'   if `stop` is provided, `duration` is recomputed.
#' @param options Named list of solver options to override defaults.
#'
#'   For ODE models (default: `list(method = "lsoda")`):
#'   - `method`: integration method passed to [deSolve::ode()]. Common
#'     choices include `"lsoda"` (default, adaptive), `"rk4"` (fixed-step
#'     Runge-Kutta), and `"euler"`.
#'
#'   For SSA models (default: `list(seed = NULL, method = "adaptive-tau")`):
#'   - `seed`: random number seed for reproducibility, or `NULL` (default)
#'     for non-deterministic results.
#'   - `method`: `"adaptive-tau"` (default, faster approximate) or `"exact"`
#'     (Gillespie's direct method).
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
#'     params = list(beta = 0.3, gamma = 0.1),
#'     equations = function() list(-beta * S * I, beta * S * I - gamma * I),
#'     times = list(duration = 50)
#' )
#' result = run_model(m)
#' result = run_model(m, params = list(beta = 0.5))
#' result = run_model(m, times = list(duration = 200))
#' result = run_model(m, options = list(method = "rk4"))
#'
#' @export
run_model = function(model, init = NULL, params = NULL, times = NULL,
    options = NULL, ...) UseMethod("run_model")

#' @export
print.modeller = function(x, ...)
{
    type = if (!is.null(x$type)) x$type else "Unknown"

    cat(type, "model\n")
    cat("  Compartments:", paste0(names(x$init), "(0) = ", x$init, collapse = ", "), "\n")
    cat("  Parameters:  ", paste0(names(x$params), " = ", x$params, collapse = ", "), "\n")
    cat("  Time:        ", x$times$start, "to", x$times$stop)
    if (!is.null(x$times$step)) cat(", step", x$times$step)
    cat("\n")
    invisible(x)
}

# Validate inputs common to all model constructors
validate_inputs = function(init, params, equations, times)
{
    if (!is.list(init) || is.null(names(init)) || any(names(init) == ""))
        stop("`init` must be a named list", call. = FALSE)
    if (!all(vapply(init, is.numeric, logical(1))))
        stop("all elements of `init` must be numeric", call. = FALSE)
    if (!is.list(params) || is.null(names(params)) || any(names(params) == ""))
        stop("`params` must be a named list", call. = FALSE)
    if (!all(vapply(params, is.numeric, logical(1))))
        stop("all elements of `params` must be numeric", call. = FALSE)
    if (!is.function(equations))
        stop("`equations` must be a function", call. = FALSE)
    if (length(formals(equations)) > 0)
        stop("`equations` must be a function with no arguments", call. = FALSE)
    if (!is.list(times) || is.null(names(init)) || any(names(init) == ""))
        stop("`times` must be a named list", call. = FALSE)
}

# Resolve times list: ensures start, stop, duration, and (for ODE) step are set.
# `overrides` is a named list of user-supplied fields that take precedence when
# determining whether stop or duration should be recomputed.
resolve_times = function(base, overrides = list(), has_step = FALSE)
{
    times = modifyList(base, overrides)
    times$start = times$start %||% 0

    if ("duration" %in% names(overrides)) {
        times$stop = times$start + times$duration
    } else if ("stop" %in% names(overrides)) {
        times$duration = times$stop - times$start
    } else if (!is.null(times$duration) && !is.null(times$stop)) {
        # Both present from base (already resolved), keep as-is
    } else if (!is.null(times$duration)) {
        times$stop = times$start + times$duration
    } else if (!is.null(times$stop)) {
        times$duration = times$stop - times$start
    } else {
        times$duration = 100
        times$stop = times$start + 100
    }

    if (has_step) times$step = times$step %||% 1
    times
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
    prefix = cumulative_prefix
    total_cols = grep(paste0("^", prefix), names(data), value = TRUE)
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
            new_name = sub(paste0("^", prefix), "", col)
            vals = summ[, col] / grid_step
            inc_data[[new_name]] = c(vals[1], diff(vals))
        }

        data[total_cols] = NULL
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
            new_name = sub(paste0("^", prefix), "", col)
            vals = data[[col]] / dt
            data[[new_name]] = c(vals[1], diff(vals))
        }
        data[total_cols] = NULL
    }

    attr(data, "dt") = dt
    attr(data, "geom") = orig_geom
    class(data) = orig_class
    data
}

# Helper to detect and issue a warning about some problem with the model
# results.
detect_problem = function(result, predicate, description, fix)
{
    problematic = apply(result[, -1], 2, predicate)
    if (any(problematic, na.rm = TRUE)) {
        first_row = which(problematic, arr.ind = TRUE)[1, 1]
        warning(description, " detected at time t = ",
            result[first_row, 1], ". ", fix, call. = FALSE)
    }
}

standard_checks = function(result)
{
    detect_problem(result, \(x) x < 0, "Negative compartment",
        "This may be because because parameters or initial conditions are negative, step size is too small, or the model logic is incorrect.")
    detect_problem(result, \(x) !is.finite(x), "Non-finite number",
        "This may be because the model is generating NAs or conducting undefined mathematical operations (such as log(-1) or division by zero).")
}
