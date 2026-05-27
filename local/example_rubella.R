# RUBELLA VACCINATION
#
# Age-structured rubella model adapted from the Berkeley Madonna model
# `rubvacc_cb.mmd` used in the modelling course (Week 2, session 10).
#
# 60 single-year age cohorts (R 1-indexed: position 1 is age 0, position 14
# is age 13). Within each year, normal SEIR-style dynamics. At year-end
# everyone ages one step; cohort 14 (age 13) receives vaccination.

library(modeller)

n_ages = 60

init = list(
    Sus    = c(1000, rep(999, n_ages - 1)),  # age 0 starts with the year's birth cohort
    Preinf = rep(0, n_ages),
    Infous = c(0, rep(1, n_ages - 1)),
    Imm    = rep(0, n_ages)
)

# vacc[k] is the proportion of the entering cohort at age k-1 that gets
# vaccinated during their move into age k-1. Only 13-year-olds vaccinated.
prop_vacc_13   = 0.7
yr_start_vacc  = 150
vacc_at_age_13 = numeric(n_ages)
vacc_at_age_13[14] = prop_vacc_13

params = list(
    R0             = 6.95,
    preinf_period  = 10,
    infous_period  = 11,
    total_popn     = 60000,
    births         = 1000,
    vacc           = vacc_at_age_13,
    yr_start_vacc  = yr_start_vacc,
    time           = c(0, 300 * 365, 1)
)

equations = function() {
    rec_rate    = 1 / infous_period
    infous_rate = 1 / preinf_period
    beta        = R0 / (total_popn * infous_period)
    foi         = beta * sum(Infous)
    year        = t / 365
    year_end    = (t %% 365) == 0
    vacc_now    = if (year > yr_start_vacc) vacc else numeric(length(vacc))

    if (year_end) {
        # Aging step. Cohort i takes the values of cohort i-1; cohort 1
        # (newborns) starts fresh from `births`. Vaccination is applied to
        # the cohort moving into age k-1: vacc_now[k] of the entering
        # cohort goes straight to Imm.
        new(Sus)    = c(births, Sus[1:59]) * (1 - vacc_now) * (1 - foi * dt)
        new(Preinf) = c(0, Preinf[1:59])
                    + foi * c(births, Sus[1:59]) * (1 - vacc_now) * dt
                    - c(0, Preinf[1:59]) * infous_rate * dt
        new(Infous) = c(0, Infous[1:59])
                    + (c(0, Preinf[1:59]) * infous_rate - c(0, Infous[1:59]) * rec_rate) * dt
        new(Imm)    = c(0, Imm[1:59])
                    + c(0, Infous[1:59]) * rec_rate * dt
                    + c(births, Sus[1:59]) * vacc_now
    } else {
        # Within-year disease dynamics, no aging
        new(Sus)    = Sus    - foi * Sus * dt
        new(Preinf) = Preinf + (foi * Sus - Preinf * infous_rate) * dt
        new(Infous) = Infous + (Preinf * infous_rate - Infous * rec_rate) * dt
        new(Imm)    = Imm    + Infous * rec_rate * dt
    }

    # Quantities of interest, in scalar form for plotting
    record(year         = year)
    record(prop_sus_5   = Sus[6]  / 1000)   # age 5  -> R index 6
    record(prop_sus_20  = Sus[21] / 1000)   # age 20 -> R index 21
    record(prop_sus_30  = Sus[31] / 1000)
    record(prop_sus_all = sum(Sus) / total_popn)
    record(infous_total = sum(Infous))
    record(new_infns_per_day = foi * sum(Sus))
}

m = difference_model(init, params, equations)

# This is a 300-year run with 1-day steps. ~3 seconds on a modern laptop.
system.time({result = run_model(m)})

show_model(m)

# show_model() doesn't yet support vector-valued compartments. Plot the
# scalar recordings directly instead.
library(ggplot2)

# Susceptibles by age, over time
ggplot(result, aes(x = year)) +
    geom_line(aes(y = prop_sus_5,   colour = "age 05")) +
    geom_line(aes(y = prop_sus_20,  colour = "age 20")) +
    geom_line(aes(y = prop_sus_30,  colour = "age 30")) +
    geom_line(aes(y = prop_sus_all, colour = "all"),  linetype = "dashed") +
    labs(y = "Proportion susceptible", colour = NULL) +
    coord_cartesian(xlim = c(0, 300))

# Total infectious per day
ggplot(result, aes(x = year, y = infous_total)) +
    geom_line() +
    labs(y = "Total infectious") +
    coord_cartesian(xlim = c(0, 300))
