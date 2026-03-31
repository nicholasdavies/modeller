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

    dS = -beta * I/N * S
    dE =  beta * I/N * S - infectious_rate * E
    dI =                   infectious_rate * E - rec_rate * I
    dR =                                         rec_rate * I
    dC =  beta * I/N * S

    return (list(dS, dE, dI, dR, dC))
}

times = list(
    start = 0,
    stop = 150
)

m = ode1_model(init, params, equations, times)
result = run_model(m)
plot(result)
show_model(m)
