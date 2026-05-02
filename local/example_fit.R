# Fitting example
library(modeller)

# Data
cases = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 4, 0, 6,
8, 15, 16, 23, 23, 29, 65, 49, 128, 154, 78, 69, 96, 102, 55,
34, 36, 19, 9, 12, 6, 5, 2, 5, 0, 5, 1, 3, 1, 1, 2, 0, 1, 0,
0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
data = data.frame(t = 0:75, cases = cases)

# Compartments with initial conditions
init = list(
    S = 99999,
    E = 0,
    I = 1,
    R = 0,
    total_infection = 0
)

# Model parameters
params = list(
    latent_period = 2,
    infectious_period = 5,
    R0 = 3, # to fit
    reporting_rate = 0.5, # to fit
    time = c(0, 75)
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
    d(total_infection) = beta * I/N * S * reporting_rate
}

# Build model
m = ode_model(init, params, equations)


# Fitting demo
objective = function(model, theta, data)
{
    params = list(
        R0 = theta[1],
        reporting_rate = theta[2]
    )

    results = run_model(m, params = params)

    dpois(data$cases, results$infection, log = TRUE)
}


fit_model(m, objective = objective, theta = c(1, 0.1), lower = c(0, 0), upper = c(10, 1), data = data)
