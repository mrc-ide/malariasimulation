
test_that('parasite drug parameters and variables correctly assigned ', {

  falc_params <- get_parameters(parasite = "falciparum")
  vivax_params <- get_parameters(parasite = "vivax")

  expect_no_error(set_drugs(falc_params, drugs = list(AL_params)))
  expect_error(set_drugs(vivax_params, drugs = list(AL_params)))

  expect_error(set_drugs(falc_params, drugs = list(CQ_PQ_params)))
  expect_no_error(set_drugs(vivax_params, drugs = list(CQ_PQ_params)))

  expect_no_error(vivax_params <- set_drugs(vivax_params, list(CQ_params, CQ_PQ_params)))

  vars <- create_variables(vivax_params)
  expect_no_error(vars["ls_drug"])
  expect_no_error(vars["ls_drug_time"])

})
