# Tests for the `time` parameter parsing inside validate_inputs().

sir_eq = function() {
    d(S) = -beta * S * I
    d(I) =  beta * S * I - gamma * I
}
sir_init = list(S = 999, I = 1)

test_that("single positive number becomes c(0, N, min(1, N))", {
    m = ode_model(sir_init, list(beta = 0.3, gamma = 0.1, time = 50), sir_eq)
    expect_equal(m$params$time, c(0, 50, 1))

    m = ode_model(sir_init, list(beta = 0.3, gamma = 0.1, time = 0.5), sir_eq)
    expect_equal(m$params$time, c(0, 0.5, 0.5))
})

test_that("two-element vector becomes c(start, stop, 1)", {
    m = ode_model(sir_init, list(beta = 0.3, gamma = 0.1, time = c(10, 50)), sir_eq)
    expect_equal(m$params$time, c(10, 50, 1))
})

test_that("two-element vector with duration < 1 clamps step to duration", {
    m = ode_model(sir_init, list(beta = 0.3, gamma = 0.1, time = c(0, 0.5)), sir_eq)
    expect_equal(m$params$time, c(0, 0.5, 0.5))
})

test_that("three-element vector clamps step to (stop - start) max", {
    m = ode_model(sir_init, list(beta = 0.3, gamma = 0.1, time = c(0, 50, 0.5)), sir_eq)
    expect_equal(m$params$time, c(0, 50, 0.5))

    m = ode_model(sir_init, list(beta = 0.3, gamma = 0.1, time = c(0, 50, 100)), sir_eq)
    expect_equal(m$params$time, c(0, 50, 50))
})

test_that("invalid time forms raise informative errors", {
    expect_error(
        ode_model(sir_init, list(beta = 0.3, gamma = 0.1, time = -5), sir_eq),
        "time"
    )
    expect_error(
        ode_model(sir_init, list(beta = 0.3, gamma = 0.1, time = c(50, 0)), sir_eq),
        "time"
    )
    expect_error(
        ode_model(sir_init, list(beta = 0.3, gamma = 0.1, time = c(0, 50, -1)), sir_eq),
        "time"
    )
    expect_error(
        ode_model(sir_init, list(beta = 0.3, gamma = 0.1, time = c(0, 1, 2, 3)), sir_eq),
        "time"
    )
})

test_that("missing time triggers warning and falls back to default", {
    expect_warning(
        m <- ode_model(sir_init, list(beta = 0.3, gamma = 0.1), sir_eq),
        "time"
    )
    expect_equal(m$params$time, c(0, 100, 1))
})
