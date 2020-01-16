test_that('discretise works on linear data', {
  continuous <- c(.1, .2, .3, .4, .5, .6, .7, .8, .9, 1.)
  expect_equal(
    discretise(continuous, 5),
    c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5)
  )
})

test_that('sample_maternal_immunity correctly samples maternal immunity from the population', {
  immunity <- c(1, 3, 7, 5, 6, 2, 3)
  sampleable <- c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE)
  groups <- c(1, 1, 1, 2, 2, 3, 3)
  died <- c(2, 4, 5)
  pm <- .5
  new_immunity <- sample_maternal_immunity(
    immunity,
    sampleable,
    died,
    groups,
    pm
  )
  expect(new_immunity[[1]] %in% c(.5, 1.5), 'first immunity valid')
  expect(new_immunity[[2]] %in% c(2.5, 3.), 'second immunity valid')
  expect(new_immunity[[3]] %in% c(2.5, 3.), 'third immunity valid')
})

test_that('sample_maternal_immunity fails if there are not enough individuals', {
  immunity <- c(1, 3, 7, 5, 6, 2, 3)
  sampleable <- c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE)
  groups <- c(1, 1, 1, 2, 2, 3, 3)
  died <- c(2, 4, 5)
  pm <- .5
  expect_error(
    sample_maternal_immunity(
      immunity,
      sampleable,
      died,
      groups,
      pm
    ),
    '*'
  )
})
