# Tests for ode_model() and run_model.ode_model().

sir_init = list(S = 999, I = 1, R = 0)
sir_params = list(beta = 0.3, gamma = 0.1, time = 50)
sir_eq = function() {
    d(S) = -beta * S * I / 1000
    d(I) =  beta * S * I / 1000 - gamma * I
    d(R) =  gamma * I
}

test_that("ode_model returns a modeller object with expected structure", {
    m = ode_model(sir_init, sir_params, sir_eq)
    expect_s3_class(m, "ode_model")
    expect_s3_class(m, "modeller")
    expect_equal(m$type, "ODE")
    expect_named(m$init, c("S", "I", "R"))
    expect_true(is.function(m$equations))
})

test_that("run_model returns a model_result data frame with t and compartments", {
    m = ode_model(sir_init, sir_params, sir_eq)
    r = run_model(m)
    expect_s3_class(r, "model_result")
    expect_s3_class(r, "data.frame")
    expect_named(r, c("t", "S", "I", "R"))
    expect_equal(r$t[1], 0)
    expect_equal(r$t[nrow(r)], 50)
    expect_equal(attr(r, "geom"), "line")
})

test_that("conservation: S + I + R is constant for closed SIR", {
    m = ode_model(sir_init, sir_params, sir_eq)
    r = run_model(m)
    totals = r$S + r$I + r$R
    expect_equal(totals, rep(1000, length(totals)), tolerance = 1e-3)
})

test_that("init/params/options can be overridden in run_model", {
    m = ode_model(sir_init, sir_params, sir_eq)
    r1 = run_model(m, params = list(beta = 0.6))
    r2 = run_model(m, params = list(beta = 0.3))
    expect_gt(max(r1$I), max(r2$I))

    r_long = run_model(m, params = list(time = 100))
    expect_equal(r_long$t[nrow(r_long)], 100)

    r_rk = run_model(m, options = list(method = "rk4"))
    expect_s3_class(r_rk, "model_result")
})

test_that("default solver method is lsode", {
    m = ode_model(sir_init, sir_params, sir_eq)
    expect_equal(m$options$method, "lsode")
})

test_that("ODE result reaches the requested stop time even when (stop-start) is not a multiple of step", {
    m = ode_model(sir_init,
        list(beta = 0.3, gamma = 0.1, time = c(0, 30.5, 1)),
        sir_eq)
    r = run_model(m)
    expect_equal(r$t[nrow(r)], 30.5)
})
