#' @export
ode_model = function(model, init, params, times)
{
    # Put model function into correct form
    # TODO check function's signature
    formals(model) = alist(t =, .state =, .params =)
    body(model) = rlang::expr({
        .x <- with(as.list(c(.state, .params)), !!body(model))
        return (list(unlist(.x)))
    })

    # Solver: returns data frame with t + compartment columns
    run_model = function(init_list, params_list, times) {
        sol = deSolve::ode(unlist(init_list), times, model, params_list)
        data = as.data.frame(sol)
        names(data)[1] = "t"
        data
    }

    model_app(run_model, init, params, times)
}
