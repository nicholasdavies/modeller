# Tests for vector-valued compartments across model types.

test_that("pack_state and unpack_state round-trip", {
    x = list(Sus = c(99, 99, 99), Inc = c(1, 0, 0), total = 0)
    flat = pack_state(x)
    expect_equal(names(flat), c("Sus.1", "Sus.2", "Sus.3", "Inc.1", "Inc.2", "Inc.3", "total"))
    schema = vapply(x, length, integer(1))
    back = unpack_state(flat, schema)
    expect_equal(back, x, ignore_attr = TRUE)
})

test_that("ODE conserves mass with vector compartments", {
    init = list(Sus = c(99, 99, 99), Inc = c(1, 0, 0), Imm = c(0, 0, 0))
    params = list(beta = 0.3, gamma = 0.1, time = c(0, 50, 1))
    eq = function() {
        N = sum(Sus + Inc + Imm)
        foi = beta * sum(Inc) / N
        d(Sus) = -foi * Sus
        d(Inc) =  foi * Sus - gamma * Inc
        d(Imm) =  gamma * Inc
    }
    m = ode_model(init, params, eq)
    r = run_model(m)
    flat_cols = setdiff(names(r), "t")
    expect_equal(flat_cols, c("Sus.1", "Sus.2", "Sus.3", "Inc.1", "Inc.2", "Inc.3", "Imm.1", "Imm.2", "Imm.3"))
    totals = rowSums(r[, flat_cols])
    expect_equal(totals, rep(298, length(totals)), tolerance = 1e-6)
})

test_that("difference equation conserves mass with vector compartments", {
    init = list(Sus = c(99, 99, 99), Inc = c(1, 0, 0), Imm = c(0, 0, 0))
    params = list(beta = 0.3, gamma = 0.1, time = c(0, 50, 1))
    eq = function() {
        N = sum(Sus + Inc + Imm)
        foi = beta * sum(Inc) / N
        new(Sus) = Sus - foi * Sus * dt
        new(Inc) = Inc + (foi * Sus - gamma * Inc) * dt
        new(Imm) = Imm + gamma * Inc * dt
    }
    m = difference_model(init, params, eq)
    r = run_model(m)
    totals = rowSums(r[, setdiff(names(r), "t")])
    expect_equal(totals, rep(298, length(totals)), tolerance = 1e-6)
})

test_that("SSA generates K transitions per vector tx() group", {
    init = list(Sus = c(99, 99, 99), Inc = c(1, 0, 0), Imm = c(0, 0, 0))
    params = list(beta = 0.3, gamma = 0.1, time = c(0, 30))
    eq = function() {
        N = sum(Sus + Inc + Imm)
        tx(Sus -> Inc) = beta * sum(Inc) / N * Sus
        tx(Inc -> Imm) = gamma * Inc
    }
    m = ssa_model(init, params, eq)
    expect_length(m$transitions, 6)   # 3 ages × 2 transition types
    # First: Sus[1] -> Inc[1]
    expect_equal(unname(m$transitions[[1]]), c(-1, 0, 0, 1, 0, 0, 0, 0, 0))
    # Fourth: Inc[1] -> Imm[1]
    expect_equal(unname(m$transitions[[4]]), c(0, 0, 0, -1, 0, 0, 1, 0, 0))
})

test_that("SSA conserves mass with vector compartments", {
    init = list(Sus = c(99, 99, 99), Inc = c(5, 0, 0), Imm = c(0, 0, 0))
    params = list(beta = 0.3, gamma = 0.1, time = c(0, 50))
    eq = function() {
        N = sum(Sus + Inc + Imm)
        tx(Sus -> Inc) = beta * sum(Inc) / N * Sus
        tx(Inc -> Imm) = gamma * Inc
    }
    m = ssa_model(init, params, eq)
    r = run_model(m, options = list(seed = 1))
    totals = rowSums(r[, setdiff(names(r), "t")])
    expect_true(all(totals == 302))
})

test_that("SSA rejects vector compartments of mismatched length", {
    init = list(Sus = c(99, 99, 99), Inc = c(1, 0))   # different lengths
    params = list(beta = 0.3, time = c(0, 30))
    eq = function() {
        tx(Sus -> Inc) = beta * Sus
    }
    expect_error(ssa_model(init, params, eq), "different.*length")
})

test_that("scalar broadcast in SSA tx() with mixed scalar+vector", {
    init = list(Sus = c(99, 99, 99), Inc = c(1, 0, 0), total_infections = 0)
    params = list(beta = 0.3, time = c(0, 30))
    eq = function() {
        # total_infections (scalar) increments on every Sus[i] -> Inc[i] event
        tx(Sus -> Inc + total_infections) = beta * Sus
    }
    m = ssa_model(init, params, eq)
    expect_length(m$transitions, 3)
    # Each transition decrements one Sus, increments one Inc, and increments total_infections
    for (i in 1:3) {
        incr = m$transitions[[i]]
        expect_equal(sum(incr[grep("^Sus", names(incr))]), -1)
        expect_equal(sum(incr[grep("^Inc", names(incr))]), 1)
        expect_equal(unname(incr["total_infections"]), 1)
    }
})

test_that("recorder works with vector compartments", {
    init = list(Sus = c(99, 99, 99), Inc = c(1, 0, 0), Imm = c(0, 0, 0))
    params = list(beta = 0.3, gamma = 0.1, time = c(0, 30, 1))
    eq = function() {
        N = sum(Sus + Inc + Imm)
        foi = beta * sum(Inc) / N
        d(Sus) = -foi * Sus
        d(Inc) =  foi * Sus - gamma * Inc
        d(Imm) =  gamma * Inc
        record(sum_inc = sum(Inc))
        record(inc_age2 = Inc[2])
    }
    m = ode_model(init, params, eq)
    r = run_model(m)
    expect_true(all(c("sum_inc", "inc_age2") %in% names(r)))
    # sum_inc equals row-wise sum of Inc.1..Inc.3
    expect_equal(r$sum_inc, r$Inc.1 + r$Inc.2 + r$Inc.3)
    # inc_age2 equals Inc.2
    expect_equal(r$inc_age2, r$Inc.2)
})

test_that("scalar models still work after vector-state refactor", {
    init = list(S = 999, I = 1, R = 0)
    params = list(beta = 0.3, gamma = 0.1, time = 30)
    eq = function() {
        d(S) = -beta * S * I / 1000
        d(I) =  beta * S * I / 1000 - gamma * I
        d(R) =  gamma * I
    }
    m = ode_model(init, params, eq)
    r = run_model(m)
    expect_equal(names(r), c("t", "S", "I", "R"))
    totals = r$S + r$I + r$R
    expect_equal(totals, rep(1000, length(totals)), tolerance = 1e-3)
})
