#' Example SIR model
#'
#' Returns an 'example' SIR ODE model for testing purposes. There is a single
#' population with homogeneous mixing, a transmission rate `beta` and a
#' recovery rate `gamma` in the model.
#'
#'
#' @param S0,I0,R0 Initial number of individuals in the S, I, and R states.
#' @param beta,gamma Transmission rate and recovery rate, respectively.
#' @param time Duration of the simulation.
#'
#' @return An object of class `c("ode_model", "modeller")`. Use
#'   [run_model()] to simulate, [plot()][plot.model_result] to visualise,
#'   and [show_model()] to launch the interactive app.
#'
#' @examples
#' m = sir_model()
#' result = run_model(m)
#' plot(result)
#'
#' @export
sir_model = function(S0 = 999, I0 = 1, R0 = 0, beta = 1, gamma = 0.5, time = 100)
{
    ode_model(
        init = list(S = S0, I = I0, R = R0),
        params = list(beta = beta, gamma = gamma, time = time),
        equations = function() {
            N = S + I + R
            d(S) = -beta * I / N * S
            d(I) =  beta * I / N * S - gamma * I
            d(R) =  gamma * I
        }
    )
}
