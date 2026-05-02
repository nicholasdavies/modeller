# Tests for total_ → incidence conversion via compute_incidence().

sir_init = list(S = 999, I = 1, R = 0, total_infection = 0)
sir_params = list(beta = 0.3, gamma = 0.1, time = c(0, 50, 1))
sir_eq = function() {
    lambda = beta * I / 1000
    d(S) = -lambda * S
    d(I) =  lambda * S - gamma * I
    d(R) =  gamma * I
    d(total_infection) = lambda * S
}

test_that("total_X column is replaced by incidence column X", {
    m = ode_model(sir_init, sir_params, sir_eq)
    r = run_model(m)
    expect_false("total_infection" %in% names(r))
    expect_true("infection" %in% names(r))
})

test_that("incidence integrates back to the cumulative total", {
    m = ode_model(sir_init, sir_params, sir_eq)
    r = run_model(m)
    # incidence is per unit time; cumulative ≈ sum * dt
    dt = attr(r, "dt")
    total_recovered = r$R[nrow(r)]
    cumulative_inc = sum(r$infection) * dt
    # Cumulative infections should be roughly R + I (everyone infected, minus initial)
    expect_equal(cumulative_inc, total_recovered + r$I[nrow(r)] - sir_init$I,
        tolerance = 1.0)
})

test_that("incidence column is non-negative", {
    m = ode_model(sir_init, sir_params, sir_eq)
    r = run_model(m)
    expect_true(all(r$infection >= -1e-6))
})
