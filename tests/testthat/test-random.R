test_that("fast_weighted_sample produces a vaguely correct distribution", {
  n <- 1e6
  prob <- c(.1, .2, .3, .4)
  counts <- as.numeric(table(fast_weighted_sample(n, prob)))
  for (i in seq_along(prob)) {
    t <- prop.test(counts[[i]], n, prob[[i]])
    actual_p <- counts[[i]] / n
    expect_gte(actual_p, t$conf.int[[1]])
    expect_lte(actual_p, t$conf.int[[2]])
  }
})

test_that("fast_weighted_sample works with un-normalised probs", {
  n <- 1e5
  prob <- c(.1, .2, .3, .4) * 100
  expected <- prob / sum(prob)
  counts <- table(fast_weighted_sample(n, prob))
  for (i in seq_along(prob)) {
    t <- prop.test(counts[[i]], n, expected[[i]])
    actual_p <- counts[[i]] / n
    expect_gte(actual_p, t$conf.int[[1]])
    expect_lte(actual_p, t$conf.int[[2]])
  }
})

test_that("fast_weighted_sample works with equal probs", {
  n <- 1e5
  prob <- c(.1, .1, .1, .1)
  expected <- prob / sum(prob)
  counts <- table(fast_weighted_sample(n, prob))
  for (i in seq_along(prob)) {
    t <- prop.test(counts[[i]], n, expected[[i]])
    actual_p <- counts[[i]] / n
    expect_gte(actual_p, t$conf.int[[1]])
    expect_lte(actual_p, t$conf.int[[2]])
  }
})

test_that("fast_weighted_sample works with light tails", {
  n <- 1e5
  prob <- c(.5, .4, .3, .2)
  expected <- prob / sum(prob)
  counts <- table(fast_weighted_sample(n, prob))
  for (i in seq_along(prob)) {
    t <- prop.test(counts[[i]], n, expected[[i]])
    actual_p <- counts[[i]] / n
    expect_gte(actual_p, t$conf.int[[1]])
    expect_lte(actual_p, t$conf.int[[2]])
  }
})
