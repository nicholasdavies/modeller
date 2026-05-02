# Tests for validate_inputs() error paths on init/params/equations shape.

ok_init = list(S = 999, I = 1)
ok_params = list(beta = 0.3, gamma = 0.1, time = 50)
ok_eq = function() {
    d(S) = -beta * S * I
    d(I) =  beta * S * I - gamma * I
}

test_that("init must be a named list of numerics", {
    expect_error(ode_model(c(S = 999, I = 1), ok_params, ok_eq), "named list")
    expect_error(ode_model(list(999, 1), ok_params, ok_eq), "named list")
    expect_error(ode_model(list(S = "x", I = 1), ok_params, ok_eq), "numeric")
})

test_that("params must be a named list of numerics", {
    expect_error(ode_model(ok_init, c(beta = 0.3), ok_eq), "named list")
    expect_error(ode_model(ok_init, list(0.3), ok_eq), "named list")
    expect_error(ode_model(ok_init, list(beta = "x", time = 10), ok_eq), "numeric")
})

test_that("equations must be a no-argument function", {
    expect_error(ode_model(ok_init, ok_params, "not a function"),
        "function")
    expect_error(ode_model(ok_init, ok_params, function(x) x),
        "no arguments")
})

test_that("ode_model rejects equations using new() or tx()", {
    bad_new = function() {
        new(S) = S - 1
        new(I) = I + 1
    }
    expect_error(ode_model(ok_init, ok_params, bad_new), "new\\(\\)")

    bad_tx = function() {
        tx(S -> I) = beta * S * I
    }
    expect_error(ode_model(ok_init, ok_params, bad_tx), "tx\\(\\)")
})

test_that("difference_model rejects equations using d() or tx()", {
    bad_d = function() {
        d(S) = -beta * S * I
        d(I) =  beta * S * I - gamma * I
    }
    expect_error(difference_model(ok_init, ok_params, bad_d), "d\\(\\)")

    bad_tx = function() {
        tx(S -> I) = beta * S * I
    }
    expect_error(difference_model(ok_init, ok_params, bad_tx), "tx\\(\\)")
})

test_that("ssa_model rejects equations using d() or new()", {
    bad_d = function() {
        d(S) = -beta * S * I
    }
    expect_error(ssa_model(ok_init, ok_params, bad_d), "d\\(\\)")

    bad_new = function() {
        new(S) = S - 1
    }
    expect_error(ssa_model(ok_init, ok_params, bad_new), "new\\(\\)")
})

test_that("ssa_model rejects malformed tx() calls", {
    # Missing the assignment
    bad_tx_no_rate = function() {
        tx(S -> I)
    }
    expect_error(ssa_model(ok_init, ok_params, bad_tx_no_rate), "tx\\(\\)")
})

test_that("init cannot use reserved or internal names", {
    expect_error(
        ode_model(list(t = 5, S = 1), ok_params, ok_eq),
        "init"
    )
    expect_error(
        ode_model(list(dt = 5, S = 1), ok_params, ok_eq),
        "init"
    )
    expect_error(
        ode_model(list(time = 5, S = 1), ok_params, ok_eq),
        "init"
    )
    expect_error(
        ode_model(list(.x = 5, S = 1), ok_params, ok_eq),
        "init"
    )
})

test_that("params cannot use reserved or internal names", {
    expect_error(
        ode_model(ok_init, list(t = 5, time = 30), ok_eq),
        "params"
    )
    expect_error(
        ode_model(ok_init, list(dt = 5, time = 30), ok_eq),
        "params"
    )
    expect_error(
        ode_model(ok_init, list(.x = 5, time = 30), ok_eq),
        "params"
    )
})

test_that("init and params cannot share names", {
    expect_error(
        ode_model(list(S = 1, beta = 999), list(beta = 0.3, time = 30), ok_eq),
        "share names"
    )
})

test_that("ode_model rejects d() referring to unknown compartments", {
    bad = function() {
        d(S) = -beta * S * I
        d(I) =  beta * S * I - gamma * I
        d(Z) = 5
    }
    expect_error(ode_model(ok_init, ok_params, bad), "Z")
})

test_that("difference_model rejects new() referring to unknown compartments", {
    bad = function() {
        new(S) = S - beta * S * I * dt
        new(I) = I + beta * S * I * dt - gamma * I * dt
        new(Z) = 0
    }
    expect_error(difference_model(ok_init, ok_params, bad), "Z")
})

test_that("ode model with missing d() leaves that compartment constant", {
    # Compartment Z is in init but has no d() — should stay constant at 5
    init = list(S = 999, I = 1, Z = 5)
    params = list(beta = 0.3, gamma = 0.1, time = 20)
    eq = function() {
        d(S) = -beta * S * I / 1000
        d(I) =  beta * S * I / 1000 - gamma * I
    }
    m = ode_model(init, params, eq)
    r = run_model(m)
    expect_true(all(abs(r$Z - 5) < 1e-8))
})

test_that("difference model with missing new() leaves that compartment constant", {
    init = list(S = 999, I = 1, Z = 7)
    params = list(beta = 0.3, gamma = 0.1, time = c(0, 20, 1))
    eq = function() {
        new(S) = S - beta * S * I * dt / 1000
        new(I) = I + beta * S * I * dt / 1000 - gamma * I * dt
    }
    m = difference_model(init, params, eq)
    r = run_model(m)
    expect_true(all(r$Z == 7))
})
