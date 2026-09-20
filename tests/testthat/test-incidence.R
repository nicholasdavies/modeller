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

test_that("results contain both total_X and its incidence column X", {
    m = ode_model(sir_init, sir_params, sir_eq)
    r = run_model(m)
    expect_named(r, c("t", "S", "I", "R", "total_infection", "infection"))
    expect_equal(r$infection, c(r$total_infection[1], diff(r$total_infection)))
})

test_that("difference and SSA results also keep total_X columns", {
    diff_eq = function() {
        lambda = beta * I / 1000
        new(S) = S - lambda * S
        new(I) = I + lambda * S - gamma * I
        new(R) = R + gamma * I
        new(total_infection) = total_infection + lambda * S
    }
    r = run_model(difference_model(sir_init, sir_params, diff_eq))
    expect_true(all(c("total_infection", "infection") %in% names(r)))

    ssa_eq = function() {
        tx(S -> I + total_infection) = beta * I / 1000 * S
        tx(I -> R) = gamma * I
    }
    r = run_model(ssa_model(sir_init, sir_params, ssa_eq),
        options = list(seed = 1))
    expect_true(all(c("total_infection", "infection") %in% names(r)))
    expect_true(all(diff(r$total_infection) >= 0))
    expect_equal(r$total_infection[nrow(r)], r$I[nrow(r)] + r$R[nrow(r)] - 1)
})

test_that("plot omits total_X columns unless requested in series", {
    m = ode_model(sir_init, sir_params, sir_eq)
    r = run_model(m)
    plotted = function(p) p$layers[[2]]$data$name

    p = plot(r)
    expect_setequal(unique(as.character(plotted(p))), c("S", "I", "R", "infection"))
    expect_equal(levels(plotted(p)), c("S", "I", "R", "infection"))

    p = plot(r, series = c("I", "total_infection"))
    expect_setequal(unique(as.character(plotted(p))), c("I", "total_infection"))
    # Requested totals go last, so other series keep their colours
    expect_equal(levels(plotted(p)), c("S", "I", "R", "infection", "total_infection"))

    p = plot(r, series = "S", vline = 10)
    expect_s3_class(p, "ggplot")
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
