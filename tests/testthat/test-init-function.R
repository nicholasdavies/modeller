# Tests for initial conditions supplied as a function of the parameters.

init_list = list(S = 990, I = 10, R = 0)
init_fn = function() list(S = N - I0, I = I0, R = 0)
params = list(N = 1000, I0 = 10, beta = 0.3, gamma = 0.1, time = 30)

ode_eq = function() {
    d(S) = -beta * I / N * S
    d(I) =  beta * I / N * S - gamma * I
    d(R) =  gamma * I
}
diff_eq = function() {
    new(S) = S - beta * I / N * S * dt
    new(I) = I + beta * I / N * S * dt - gamma * I * dt
    new(R) = R + gamma * I * dt
}
ssa_eq = function() {
    tx(S -> I) = beta * I / N * S
    tx(I -> R) = gamma * I
}

test_that("list and function forms of init give identical results", {
    expect_equal(run_model(ode_model(init_fn, params, ode_eq)),
        run_model(ode_model(init_list, params, ode_eq)))
    expect_equal(run_model(difference_model(init_fn, params, diff_eq)),
        run_model(difference_model(init_list, params, diff_eq)))
    expect_equal(
        run_model(ssa_model(init_fn, params, ssa_eq), options = list(seed = 1)),
        run_model(ssa_model(init_list, params, ssa_eq), options = list(seed = 1)))
})

test_that("init function fixes the model structure at default parameters", {
    m = ode_model(function() list(S = N - I0, I = I0, R = 0, total_inf = 0),
        params, function() {
            d(S) = -beta * I / N * S
            d(I) =  beta * I / N * S - gamma * I
            d(R) =  gamma * I
            d(total_inf) = beta * I / N * S
        })
    expect_equal(m$init, list(S = 990, I = 10, R = 0, total_inf = 0))
    expect_true(is.function(m$init_fn))
    expect_true("inf" %in% names(run_model(m)))
    expect_null(ode_model(init_list, params, ode_eq)$init_fn)
})

test_that("init function is re-evaluated with overridden parameters", {
    m = ode_model(init_fn, params, ode_eq)
    r = run_model(m, params = list(I0 = 100, N = 2000))
    expect_equal(r$S[1], 1900)
    expect_equal(r$I[1], 100)
})

test_that("init function sees t as the start time", {
    m = ode_model(function() list(S = 1000 - t, I = t, R = 0),
        modifyList(params, list(time = c(5, 30))), ode_eq)
    expect_equal(m$init, list(S = 995, I = 5, R = 0))
    r = run_model(m, params = list(time = c(8, 30)))
    expect_equal(r$I[1], 8)
})

test_that("init function can use variables from its enclosing environment", {
    scale = 2
    m = ode_model(function() list(S = N - I0 * scale, I = I0 * scale, R = 0),
        params, ode_eq)
    expect_equal(m$init$I, 20)
})

test_that("init function may reorder compartments", {
    m = ode_model(init_fn, params, ode_eq)
    m$init_fn = build_init_fn(function() list(R = 0, I = I0, S = N - I0))
    r = run_model(m)
    expect_equal(r$S[1], 990)
    expect_equal(r$I[1], 10)
})

test_that("a list passed to run_model overrides the init function", {
    m = ode_model(init_fn, params, ode_eq)
    r = run_model(m, init = list(I = 50), params = list(I0 = 20))
    expect_equal(r$S[1], 980)
    expect_equal(r$I[1], 50)
})

test_that("a function passed to run_model overrides a subset of compartments", {
    m = ode_model(init_list, params, ode_eq)
    r = run_model(m, init = function() list(I = I0 * 3), params = list(I0 = 4))
    expect_equal(r$S[1], 990)
    expect_equal(r$I[1], 12)
})

test_that("run_model rejects init overrides incompatible with the model", {
    m = ode_model(init_fn, params, ode_eq)
    expect_error(run_model(m, init = function() list(Z = 1)), "unknown compartment")
    expect_error(run_model(m, init = list(Z = 1)), "unknown compartment")
    expect_error(run_model(m, init = function() list(I = c(1, 2))), "wrong length")
    expect_error(run_model(m, init = function() c(I = 1)), "named list")
    expect_error(run_model(m, init = function(x) list(I = 1)), "no arguments")
})

test_that("model's init function must keep returning the same compartments", {
    m = ode_model(function() if (I0 > 50) list(S = N - I0, I = I0) else
            list(S = N - I0, I = I0, R = 0),
        params, ode_eq)
    expect_error(run_model(m, params = list(I0 = 100)), "did not return")
})

test_that("constructors validate init functions", {
    expect_error(ode_model(function(x) list(S = 1), params, ode_eq), "no arguments")
    expect_error(ode_model(function() c(S = 1), params, ode_eq), "named list")
    expect_error(ode_model(function() list(S = "a"), params, ode_eq), "numeric")
    expect_error(ode_model(function() list(S = 1, beta = 2), params, ode_eq),
        "share names")
})

test_that("print notes an init function", {
    expect_output(print(ode_model(init_fn, params, ode_eq)), "from init function")
    expect_no_match(capture.output(print(ode_model(init_list, params, ode_eq))),
        "from init function")
})

test_that("SSA rejects non-integer initial values from an init function", {
    m = ssa_model(function() list(S = N * 0.995, I = N * 0.005, R = 0),
        params, ssa_eq)
    expect_error(run_model(m, params = list(N = 1001), options = list(seed = 1)),
        "must be an integer")
    expect_error(run_model(m, params = list(N = 1001),
        options = list(seed = 1, method = "exact")), "must be an integer")
    expect_no_error(run_model(m, options = list(seed = 1)))
})

test_that("fit_model can estimate a parameter that sets initial conditions", {
    m = ode_model(init_fn, params, ode_eq)
    truth = run_model(m, params = list(I0 = 25))
    objective = function(model, theta, data) {
        results = run_model(model, params = list(I0 = theta[1]))
        -sum((results$I - data$I)^2)
    }
    fit = fit_model(m, objective = objective,
        theta = 10, lower = 1, upper = 100, data = truth, maxit = 200)
    expect_equal(fit$theta, 25, tolerance = 0.01)
})
