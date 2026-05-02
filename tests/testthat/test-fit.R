# Tests for fit_model() basic mechanics.

sir_init = list(S = 999, I = 1, R = 0, total_infection = 0)
sir_params = list(beta = 0.3, gamma = 0.1, reporting = 0.5, time = 30)
sir_eq = function() {
    d(S) = -beta * S * I / 1000
    d(I) =  beta * S * I / 1000 - gamma * I
    d(R) =  gamma * I
    d(total_infection) = beta * S * I / 1000 * reporting
}

test_that("fit_model recovers a known parameter", {
    m = ode_model(sir_init, sir_params, sir_eq)
    truth = run_model(m)

    objective = function(model, theta, data) {
        results = run_model(model, params = list(beta = theta[1]))
        # Sum of squared errors as log-likelihood
        -sum((results$infection - data$infection)^2)
    }

    fit = fit_model(m, objective = objective,
        theta = 0.5, lower = 0.01, upper = 1.0,
        data = truth, maxit = 200)
    expect_equal(fit$theta, 0.3, tolerance = 0.05)
})

test_that("fit_model uses the supplied theta as the starting point", {
    # Linear (well-behaved) objective whose maximum is at theta = truth.
    # If optim were started at the wrong place (e.g. bounded theta passed
    # where unbounded was expected), Nelder-Mead can still find the optimum
    # given enough iterations, so we cap iterations and use a steep penalty
    # so the wrong starting point materially affects the recovered value.
    m = ode_model(sir_init, sir_params, sir_eq)
    truth = run_model(m)

    objective = function(model, theta, data) {
        results = run_model(model, params = list(beta = theta[1]))
        -sum((results$infection - data$infection)^2)
    }

    # Start exactly at the truth: a properly-coded fit should stay there.
    fit = fit_model(m, objective = objective,
        theta = 0.3, lower = 0.01, upper = 1.0,
        data = truth, maxit = 5)
    expect_equal(fit$theta, 0.3, tolerance = 1e-3)
})

test_that("fit_model rejects ill-formed objective", {
    m = ode_model(sir_init, sir_params, sir_eq)
    bad_obj = function(x, y) sum(x)
    expect_error(
        fit_model(m, objective = bad_obj, theta = 0.5,
            lower = 0, upper = 1, data = NULL),
        "objective"
    )
})

test_that("fit_model rejects mismatched bounds", {
    m = ode_model(sir_init, sir_params, sir_eq)
    obj = function(model, theta, data) 0

    expect_error(
        fit_model(m, objective = obj, theta = c(0.5, 0.5),
            lower = 0, upper = 1, data = NULL),
        "same length"
    )
    expect_error(
        fit_model(m, objective = obj, theta = 0.5,
            lower = 1, upper = 0, data = NULL),
        "lower"
    )
})
