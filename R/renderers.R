create_renderers <- function(individuals, states, variables, parameters) {
  list(
    function(api) {
      source_humans <- api$get_state(
        individuals$human,
        states$S,
        states$U,
        states$A,
        states$D
      )
      source_mosquitos <- api$get_state(individuals$mosquito, states$Im)
      age <- api$get_variable(individuals$human, variables$age)
      ib <- api$get_variable(individuals$human, variables$ib)[source_humans]

      epsilon <- eir(
        age[source_humans],
        api$get_variable(individuals$human, variables$xi)[source_humans],
        api$get_variable(
          individuals$mosquito,
          variables$mosquito_variety
        )[source_mosquitos],
        api$get_parameters()
      )
      list(EIR = mean(epsilon))
    },
    function(api) {
      parameters <- api$get_parameters()
      source_mosquitos <- api$get_state(individuals$mosquito, states$Sm)

      age <- api$get_variable(individuals$human, variables$age)
      asymptomatic <- api$get_state(individuals$human, states$A)

      a_infectivity <- asymptomatic_infectivity(
        age[asymptomatic],
        api$get_variable(individuals$human, variables$id)[asymptomatic],
        parameters
      )

      # Create a dataframe frame with human age, xi and infectivity
      infectivity_frame <- create_infectivity_frame(
        age,
        api$get_variable(individuals$human, variables$xi),
        list(
          list(api$get_state(individuals$human, states$D), parameters$cd),
          list(asymptomatic, a_infectivity),
          list(api$get_state(individuals$human, states$U), parameters$cu)
        )
      )

      lambda <- mosquito_force_of_infection(
        api$get_variable(
          individuals$mosquito,
          variables$mosquito_variety
        )[source_mosquitos],
        infectivity_frame,
        parameters
      )
      list(FOIM = mean(lambda))
    }
  )
}
