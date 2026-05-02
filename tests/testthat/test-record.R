# Tests for the record() feature inside equations.

sir_init = list(S = 999, I = 1, R = 0)
sir_params = list(beta = 0.3, gamma = 0.1, time = 30)

test_that("record() values appear as columns in ODE result", {
    eq = function() {
        d(S) = -beta * S * I / 1000
        d(I) =  beta * S * I / 1000 - gamma * I
        d(R) =  gamma * I
        record(half_R = R / 2)
        record(force = beta * I / 1000)
    }
    m = ode_model(sir_init, sir_params, eq)
    r = run_model(m)
    expect_true(all(c("half_R", "force") %in% names(r)))
    expect_equal(r$half_R, r$R / 2)
})

test_that("record() values respect time-dependent code paths", {
    eq = function() {
        if (t <= 10) b = beta else b = beta * 0.5
        d(S) = -b * S * I / 1000
        d(I) =  b * S * I / 1000 - gamma * I
        d(R) =  gamma * I
        record(beta_used = b)
    }
    m = ode_model(sir_init, sir_params, eq)
    r = run_model(m)
    expect_equal(r$beta_used[r$t <= 10], rep(0.3, sum(r$t <= 10)))
    expect_equal(r$beta_used[r$t > 10], rep(0.15, sum(r$t > 10)))
})

test_that("record() works for SSA models", {
    tx_eq = function() {
        tx(S -> I) = beta * S * I / 1000
        tx(I -> R) = gamma * I
        record(half_R = R / 2)
    }
    m = ssa_model(sir_init, sir_params, tx_eq)
    r = run_model(m, options = list(seed = 1))
    expect_true("half_R" %in% names(r))
    expect_equal(r$half_R, r$R / 2)
    expect_equal(attr(r, "geom"), "step")
})

test_that("record() in SSA gets the absolute t (not double-offset)", {
    tx_eq = function() {
        tx(S -> I) = beta * S * I / 1000
        tx(I -> R) = gamma * I
        record(t_recorded = t)
    }
    init = list(S = 999, I = 1, R = 0)
    p = list(beta = 0.3, gamma = 0.1, time = c(10, 50))  # start at t = 10
    m = ssa_model(init, p, tx_eq)
    r = run_model(m, options = list(seed = 1))
    expect_equal(r$t_recorded, r$t)
})

test_that("record() works for difference models", {
    diff_eq = function() {
        lambda = beta * I / 1000
        new(S) = S - lambda * S * dt
        new(I) = I + lambda * S * dt - gamma * I * dt
        new(R) = R + gamma * I * dt
        record(half_R = R / 2)
    }
    m = difference_model(sir_init, sir_params, diff_eq)
    r = run_model(m)
    expect_true("half_R" %in% names(r))
    expect_equal(r$half_R, r$R / 2)
})

test_that("record() in difference models has access to dt", {
    diff_eq = function() {
        lambda = beta * I / 1000
        new(S) = S - lambda * S * dt
        new(I) = I + lambda * S * dt - gamma * I * dt
        new(R) = R + gamma * I * dt
        record(step = dt)
    }
    p = list(beta = 0.3, gamma = 0.1, time = c(0, 30, 0.5))
    m = difference_model(sir_init, p, diff_eq)
    r = run_model(m)
    expect_equal(r$step, rep(0.5, nrow(r)))
})

test_that("model without record() has no recorder and no extra columns", {
    eq = function() {
        d(S) = -beta * S * I / 1000
        d(I) =  beta * S * I / 1000 - gamma * I
        d(R) =  gamma * I
    }
    m = ode_model(sir_init, sir_params, eq)
    expect_null(m$recorder)
    r = run_model(m)
    expect_named(r, c("t", "S", "I", "R"))
})
