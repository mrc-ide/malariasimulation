test_that("History in R", {
  h <- create_history(size = 2, 10000)
  history_push(h, 3000, 100)
  history_push(h, 5000, 200)
  
  # With linear = FALSE
  expect_equal(history_at(h, 1, FALSE), 10000) # Before first user input: Expect default 10000
  expect_equal(history_at(h, 100, FALSE), 3000) # At first user input: Expect user input
  expect_equal(history_at(h, 150, FALSE), 3000) # Between first and second user input: Expect Expect user input
  expect_equal(history_at(h, 200, FALSE), 5000) # At second user input: Expect Expect user input
  expect_equal(history_at(h, 300, FALSE), 5000) # After second user input: Expect Expect user input
  
  # With linear = TRUE
  expect_equal(history_at(h, 1, TRUE), 10000) # Before first user input: Expect default 10000
  expect_equal(history_at(h, 100, TRUE), 3000) # At first user input: Expect Expect user input
  expect_equal(history_at(h, 150, TRUE), 4000) # Between first and second user input: Expect Interpolation
  expect_equal(history_at(h, 200, TRUE), 5000) # At second user input: Expect Expect user input
  expect_equal(history_at(h, 300, TRUE), 10000) # After second user input: Expect default
})
