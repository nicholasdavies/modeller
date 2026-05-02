# Tests for ssa_model() and run_model.ssa_model().

sir_init = list(S = 999, I = 1, R = 0)
sir_params = list(beta = 0.3, gamma = 0.1, time = 50)
sir_tx = function() {
    tx(S -> I) = beta * S * I / 1000
    tx(I -> R) = gamma * I
}

test_that("ssa_model returns a modeller object with parsed transitions", {
    m = ssa_model(sir_init, sir_params, sir_tx)
    expect_s3_class(m, "ssa_model")
    expect_s3_class(m, "modeller")
    expect_equal(m$type, "SSA")
    expect_length(m$transitions, 2)
    expect_equal(m$transitions[[1]][["S"]], -1)
    expect_equal(m$transitions[[1]][["I"]], +1)
    expect_equal(m$transitions[[2]][["I"]], -1)
    expect_equal(m$transitions[[2]][["R"]], +1)
})

test_that("default options include epsilon = 0.05", {
    m = ssa_model(sir_init, sir_params, sir_tx)
    expect_equal(m$options$epsilon, 0.05)
    expect_equal(m$options$method, "adaptive-tau")
})

test_that("run_model returns step-geom result with conserved population", {
    m = ssa_model(sir_init, sir_params, sir_tx)
    r = run_model(m, options = list(seed = 1))
    expect_s3_class(r, "model_result")
    expect_equal(attr(r, "geom"), "step")
    expect_equal(attr(r, "dt"), 0)
    expect_true(all(r$S + r$I + r$R == 1000))
})

test_that("seed produces reproducible results", {
    m = ssa_model(sir_init, sir_params, sir_tx)
    r1 = run_model(m, options = list(seed = 42))
    r2 = run_model(m, options = list(seed = 42))
    expect_equal(r1$t, r2$t)
    expect_equal(r1$I, r2$I)
})

test_that("malformed transitions raise informative errors", {
    bad_tx = function() {
        tx(S + 5 -> I) = beta * S * I
    }
    expect_error(
        ssa_model(sir_init, sir_params, bad_tx),
        "transition"
    )

    unknown_tx = function() {
        tx(S -> X) = beta * S
    }
    expect_error(
        ssa_model(sir_init, sir_params, unknown_tx),
        "unknown compartment"
    )
})
