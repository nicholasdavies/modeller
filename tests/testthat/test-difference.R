# Tests for difference_model() and run_model.difference_model().

sir_init = list(S = 999, I = 1, R = 0)
sir_params = list(beta = 0.3, gamma = 0.1, time = c(0, 50, 1))
sir_eq = function() {
    lambda = beta * I / 1000
    new(S) = S - lambda * S * dt
    new(I) = I + lambda * S * dt - gamma * I * dt
    new(R) = R + gamma * I * dt
}

test_that("difference_model returns a modeller object", {
    m = difference_model(sir_init, sir_params, sir_eq)
    expect_s3_class(m, "difference_model")
    expect_s3_class(m, "modeller")
    expect_equal(m$type, "difference")
})

test_that("run_model produces step-geom result with t at each step", {
    m = difference_model(sir_init, sir_params, sir_eq)
    r = run_model(m)
    expect_s3_class(r, "model_result")
    expect_equal(attr(r, "geom"), "step")
    expect_equal(attr(r, "dt"), 1)
    expect_equal(r$t, seq(0, 50, 1))
})

test_that("conservation: total population stays at 1000", {
    m = difference_model(sir_init, sir_params, sir_eq)
    r = run_model(m)
    totals = r$S + r$I + r$R
    expect_equal(totals, rep(1000, length(totals)), tolerance = 1e-6)
})

test_that("smaller dt and explicit dt scaling give the same trajectory", {
    m1 = difference_model(sir_init,
        list(beta = 0.3, gamma = 0.1, time = c(0, 30, 1)), sir_eq)
    m2 = difference_model(sir_init,
        list(beta = 0.3, gamma = 0.1, time = c(0, 30, 0.1)), sir_eq)
    r1 = run_model(m1)
    r2 = run_model(m2)
    # At common time points, results should be similar (small Euler error)
    common = which(r2$t %in% r1$t)
    expect_equal(r2$I[common], r1$I, tolerance = 0.5)
})
