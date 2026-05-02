# Helper to transform theta from its natural bounded scale (each element
# between corresponding entries in (lower, upper) to an unbounded scale for
# fitting.
to_unbounded = function(theta, lower, upper)
{
    for (i in seq_along(theta)) {
        if (lower[[i]] != -Inf && upper[[i]] != Inf) {
            # Two bounds
            theta[[i]] = qlogis((theta[[i]] - lower[[i]]) / (upper[[i]] - lower[[i]]))
        } else if (lower[[i]] == -Inf && upper[[i]] != Inf) {
            # Upper bound only
            theta[[i]] = log(upper[[i]] - theta[[i]])
        } else if (lower[[i]] != -Inf && upper[[i]] == Inf) {
            # Lower bound only
            theta[[i]] = log(theta[[i]] - lower[[i]])
        }
    }
    return (theta)
}

# Helper to reverse the to_unbounded transformation.
to_bounded = function(theta, lower, upper)
{
    for (i in seq_along(theta)) {
        if (lower[[i]] != -Inf && upper[[i]] != Inf) {
            # Two bounds
            theta[[i]] = lower[[i]] + (upper[[i]] - lower[[i]]) * plogis(theta[[i]])
        } else if (lower[[i]] == -Inf && upper[[i]] != Inf) {
            # Upper bound only
            theta[[i]] = upper[[i]] - exp(theta[[i]])
        } else if (lower[[i]] != -Inf && upper[[i]] == Inf) {
            # Lower bound only
            theta[[i]] = exp(theta[[i]]) + lower[[i]]
        }
    }
    return (theta)
}

# Wrapper for the user's provided objective function compatible with optim.
modeller_objective = function(theta, user_model, user_data, user_objective,
    user_lower, user_upper)
{
    # Convert theta to natural scale.
    theta_bounded = to_bounded(theta, user_lower, user_upper)

    # Sum of values returned by user_objective (to maximize). This is most
    # appropriate if user_objective returns (pointwise) log-likelihoods.
    sum(user_objective(user_model, theta_bounded, user_data))
}

#' Fit a model to data
#'
#' Fits model parameters by minimising the negative sum of values returned
#' by `objective`, using Nelder-Mead optimisation via [optim()]. Parameters
#' are internally transformed to an unbounded scale so that box constraints
#' are respected without requiring a bounded optimiser.
#'
#' @param model A model object created by [ode_model()] or [ssa_model()].
#' @param objective A function with signature `function(model, theta, data)`
#'   that returns a numeric vector of pointwise log-likelihoods (or any
#'   quantity to be summed and maximised).
#' @param theta Numeric vector of initial parameter guesses. Must lie
#'   strictly between `lower` and `upper` (not on the bounds).
#' @param lower Numeric vector of lower bounds, the same length as `theta`.
#'   Use `-Inf` for unbounded below.
#' @param upper Numeric vector of upper bounds, the same length as `theta`.
#'   Use `Inf` for unbounded above.
#' @param data Data to pass to `objective` (typically a data frame of
#'   observations).
#' @param maxit Maximum number of iterations (default 500).
#'
#' @return A list with components:
#'   \describe{
#'     \item{theta}{Fitted parameter values on the original scale.}
#'     \item{value}{Maximised objective (sum of values returned by
#'       `objective`).}
#'     \item{iterations}{Number of function evaluations.}
#'     \item{convergence}{Character string describing convergence status.}
#'   }
#'
#' @export
fit_model = function(model, objective, theta, lower, upper, data, maxit = 500)
{
    # Check parameters
    if (!is.function(objective) || !identical(names(formals(objective)), c("model", "theta", "data"))) {
        stop("`objective` must be a function with parameters `model, theta, data`.")
    }
    if (!is.numeric(theta) || length(theta) == 0) {
        stop("`theta` must be a numeric vector of length one or greater.")
    }
    if (!is.numeric(upper) || !is.numeric(lower) ||
            length(lower) != length(theta) || length(upper) != length(theta)) {
        stop("`upper` and `lower` must be numeric vectors of the same length as `theta`.")
    }
    if (any(lower >= upper)) {
        stop("All values of `lower` must be lower than the corresponding value of `upper`.")
    }

    # Put initial guess for theta on the unbounded scale.
    theta_unbounded = to_unbounded(theta, lower, upper)

    # Check for any incompatibilities between theta and the bounds.
    if (any(!is.finite(theta_unbounded))) {
        stop("Initial guess for `theta` cannot lie exactly on the bounds or outside the bounds.")
    }

    # Fit the model.
    fit = optim(theta_unbounded, fn = modeller_objective, user_model = model,
        user_data = data, user_objective = objective,
        user_lower = lower, user_upper = upper,
        control = list(warn.1d.NelderMead = FALSE,
            maxit = maxit, fnscale = -1))

    # Process convergence results and any message.
    conv = "Converged."
    if (fit$convergence == 1) {
        conv = "Reached maximum number of iterations."
    } else if (fit$convergence == 10) {
        conv = "Nelder-Mead simplex degenerated."
    }
    if (!is.null(fit$message)) {
        conv = paste(conv, "Message:", fit$message)
    }

    return (list(theta = to_bounded(fit$par, lower, upper), value = fit$value,
        iterations = fit$counts[["function"]], convergence = conv))
}
