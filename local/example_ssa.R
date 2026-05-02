# SSA EXAMPLE
# X click on plot gives total_infections, not infections
# X click line at far left isn't visible
# X click line on step models -- not clear where measuring
# X look for d(), tx(), new() and correct user if they're using the wrong type.
# X Correct erroneous uses of d(), tx(), new().
# X Check for name clashes in init and params. And for t, dt.
# X make sure no tx/d/new are missing.

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
    R0 = 13,
    time = 150
)

equations = function()
{
    infectious_rate = 1 / latent_period
    rec_rate = 1 / infectious_period
    beta = R0 / infectious_period
    N = S + E + I + R

    tx(S -> E + total_infections) = beta * I/N * S
    tx(E -> I) = infectious_rate * E
    tx(I -> R) = rec_rate * I

    record(something = R/2)
}

m = ssa_model(init, params, equations)
result = run_model(m)
plot(result)
show_model(m)

