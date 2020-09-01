test_that('discretise normal works outside the boundaries', {
  expect_equal(
    discretise_normal(c(-100, 2, 100), 5),
    c(1, 4, 5)
  )
})
