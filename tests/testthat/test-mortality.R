test_that('sample_mothers correctly samples mothers from the population', {
  sampleable <- c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE)
  groups <- c(1, 1, 1, 2, 2, 3, 3)
  died <- c(2, 4, 5)
  mothers <- sample_mothers(
    sampleable,
    died,
    groups
  )
  expect(mothers[[1]] %in% c(1, 2), 'first immunity valid')
  expect(mothers[[2]] %in% c(4, 5), 'second immunity valid')
  expect(mothers[[3]] %in% c(4, 5), 'third immunity valid')
})

test_that('sample_mothers fails if there are not enough individuals', {
  sampleable <- c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE)
  groups <- c(1, 1, 1, 2, 2, 3, 3)
  died <- c(2, 4, 5)
  expect_error(
    sample_mothers(
      sampleable,
      died,
      groups
    ),
    '*'
  )
})
