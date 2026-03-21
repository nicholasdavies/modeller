# SSA EXAMPLE

library(modeller)

init = list(
    S = 999,
    E = 0,
    I = 1,
    R = 0,
    total_infections = 0
)

params = list(
    latent_period = 8,
    infectious_period = 7,
    R0 = 13
)

transitions = list(
    c(S = -1, E = +1, total_infections = +1),
    c(E = -1, I = +1),
    c(I = -1, R = +1)
)

model = function()
{
    infectious_rate = 1 / latent_period
    rec_rate = 1 / infectious_period
    beta = R0 / infectious_period
    N = S + E + I + R

    xS_E = beta * I/N * S
    xE_I = infectious_rate * E
    xI_R = rec_rate * I

    return (list(xS_E, xE_I, xI_R))
}

times = list(
    duration = 150
)

ssa_model(init, params, transitions, model, times)

