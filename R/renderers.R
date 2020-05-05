create_renderers <- function(individuals, states, variables, parameters) {
  list(
    function(api) {
      list(
        mean_icm=mean(api$get_variable(individuals$human, variables$icm)),
        mean_ivm=mean(api$get_variable(individuals$human, variables$ivm)),
        mean_ib=mean(api$get_variable(individuals$human, variables$ib)),
        mean_ica=mean(api$get_variable(individuals$human, variables$ica)),
        mean_iva=mean(api$get_variable(individuals$human, variables$iva))
      )
    },
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
    }
  )
}
