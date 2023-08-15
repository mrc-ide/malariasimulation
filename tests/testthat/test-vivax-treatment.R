test_that('vivax drug parameters and variables can be initiated correctly', {

  parameters <- get_parameters(parasite = "vivax")
  expect_error(set_drugs(parameters, list(CQ_params, DHA_PQP_params)))
  expect_no_error(set_drugs(parameters, list(CQ_params)))
  expect_no_error(set_drugs(parameters, list(CQ_params)))
  expect_no_error(parameters <- set_drugs(parameters, list(CQ_params, CQ_PQ_params)))

  ## Do my drug variables get made?
  vars <- create_variables(parameters)
  expect_no_error(vars["ls_drug"])
  expect_no_error(vars["ls_drug_time"])

})
