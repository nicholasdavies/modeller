library(modeller)

init = list(
    S = 99999,
    E = 0,
    I = 1,
    R = 0
)

params = list(
    latent_period = 8,
    infectious_period = 7,
    R0 = 13
)

model = function()
{
    infectious_rate = 1 / latent_period
    rec_rate = 1 / infectious_period
    beta = R0 / infectious_period
    N = S + E + I + R

    dS = -beta * I/N * S
    dE =  beta * I/N * S - E * infectious_rate
    dI =                   E * infectious_rate - I * rec_rate
    dR =                                         I * rec_rate

    return (list(dS, dE, dI, dR))
}

times = list(
    start = 0,
    stop = 150
)

ode_model(model, init, params, times)

