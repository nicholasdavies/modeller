# X multiple plots? -- in code?
# X fitting -- in code
# X euler timestep blows up axes with large time step
# X don't hide negative outputs.
# X Detect unusual outputs
# X Handle warnings or errors
# X Conversion of total_xxx should be part of running the model, not just plotting it
# X record user values within equations -- add to plot in single column
# X compare button, like in covidm_shiny2 app.
# X axis titles, plot title
# use arrows
# contact volunteers

# SEIR model for measles
library(modeller)

# Compartments with initial conditions
init <- list(
    S = 99999,
    E = 0,
    I = 1,
    R = 0,
    total_infection = 0
)

# Model parameters
params <- list(
    latent_period = 8,
    infectious_period = 7,
    R0 = 13,
    time = 150
)

# Model equations
equations <- function()
{
    infectious_rate <- 1 / latent_period
    rec_rate <- 1 / infectious_period
    beta <- R0 / infectious_period
    N <- S + E + I + R

    d(S) <- -beta * I/N * S
    d(E) <-  beta * I/N * S - infectious_rate * E
    d(I) <-                   infectious_rate * E - rec_rate * I
    d(R) <-                                         rec_rate * I
    d(total_infection) <- beta * I/N * S

    record(R0 = R0)
}

# Build model
m <- ode_model(init, params, equations)

## Run and plot model
#result <- run_model(m)
#library(ggplot2)
#plot(result) + labs(y = "Number of people")
## plot recorded variables in grid?
#cowplot::plot_grid(plot(result), plot(result), ncol = 1)

# Model explorer
show_model(m)
