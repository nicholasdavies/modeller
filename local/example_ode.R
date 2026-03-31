# SEIR model for measles
library(modeller)

# Compartments with initial conditions
init = list(
    S = 99999,
    E = 0,
    I = 1,
    R = 0
)

# Model parameters
params = list(
    latent_period = 8,
    infectious_period = 7,
    R0 = 13
)

# Model equations
equations = function()
{
    infectious_rate = 1 / latent_period
    rec_rate = 1 / infectious_period
    beta = R0 / infectious_period
    N = S + E + I + R

    d(S) = -beta * I/N * S
    d(E) =  beta * I/N * S - infectious_rate * E
    d(I) =                   infectious_rate * E - rec_rate * I
    d(R) =                                         rec_rate * I
}

# Times
times = list(
    start = 0,
    stop = 150
)

# Build model
m = ode_model(init, params, equations, times)

# Run and plot model
result = run_model(m)
plot(result)

# Model explorer
show_model(m)
