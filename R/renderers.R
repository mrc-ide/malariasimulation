create_renderers <- function(individuals, states, variables, parameters) {
  list(
    function(frame) {
      source_humans <- frame$get_state(
        individuals$human,
        states$S,
        states$U,
        states$A,
        states$D
      )
      source_mosquitos <- frame$get_state(individuals$mosquito, states$Im)
      age <- frame$get_variable(individuals$human, variables$age)
      ib <- frame$get_variable(individuals$human, variables$ib)[source_humans]

      epsilon <- eir(
        age[source_humans],
        frame$get_variable(individuals$human, variables$xi)[source_humans],
        frame$get_variable(
          individuals$mosquito,
          variables$mosquito_variety
        )[source_mosquitos],
        parameters
      )
      list(EIR = mean(epsilon))
    }
  )
}
