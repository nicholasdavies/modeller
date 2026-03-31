# SSA EXAMPLE
# TODO click on plot gives total_infections, not infections
# TODO click line at far left isn't visible
# TODO click line on step models -- not clear where measuring
# TODO look for d(), tx(), new() and correct user if they're using the wrong
#      type.
# TODO Correct erroneous uses of d(), tx(), new().
# TODO Check for name clashes in init and params. And for t, dt.
# TODO make sure no tx/d/new are missing.
# TODO everything TODO in ssa.R

library(modeller)

init = list(
    S = 999,
    E = 0,
    I = 1,
    R = 0
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

    tx(S -> E) = beta * I/N * S # TODO total_infections
    tx(E -> I) = infectious_rate * E
    tx(I -> R) = rec_rate * I
}

times = list(
    duration = 150
)

m = ssa_model(init, params, equations, times)
result = run_model(m)
plot(result)
show_model(m)

