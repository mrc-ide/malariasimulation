test_that("fast_weighted_sample produces a vaguely correct distribution", {
  n <- 1e5
  prob <- c(.1, .2, .3, .4)
  counts <- table(fast_weighted_sample(n, prob))
  for (i in seq_along(prob)) {
    t <- prop.test(counts[[i]], n)
    expect_gte(prob[[i]], t$conf.int[[1]])
    expect_lte(prob[[i]], t$conf.int[[2]])
  }
})

test_that("fast_weighted_sample works with un-normalised probs", {
  n <- 1e5
  prob <- c(.1, .2, .3, .4) * 100
  expected <- prob / sum(prob)
  counts <- table(fast_weighted_sample(n, prob))
  for (i in seq_along(prob)) {
    t <- prop.test(counts[[i]], n)
    expect_gte(expected[[i]], t$conf.int[[1]])
    expect_lte(expected[[i]], t$conf.int[[2]])
  }
})

test_that("fast_weighted_sample works with equal probs", {
  n <- 1e5
  prob <- c(.1, .1, .1, .1)
  expected <- prob / sum(prob)
  counts <- table(fast_weighted_sample(n, prob))
  for (i in seq_along(prob)) {
    t <- prop.test(counts[[i]], n)
    expect_gte(expected[[i]], t$conf.int[[1]])
    expect_lte(expected[[i]], t$conf.int[[2]])
  }
})
