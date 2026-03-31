# DIFFERENCE EQUATION EXAMPLE

library(modeller)

init = list(
    S = 99999,
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

equations = function()
{
    infectious_rate = 1 / latent_period
    rec_rate = 1 / infectious_period
    beta = R0 / infectious_period
    N = S + E + I + R

    infections = beta * I/N * S * dt
    becoming_infectious = infectious_rate * E * dt
    recoveries = rec_rate * I * dt

    list(
        S = S - infections,
        E = E + infections - becoming_infectious,
        I = I + becoming_infectious - recoveries,
        R = R + recoveries,
        total_infections = total_infections + infections
    )
}

times = list(
    duration = 150,
    step = 1
)

m = difference_model(init, params, equations, times)
result = run_model(m)
plot(result)
show_model(m)
