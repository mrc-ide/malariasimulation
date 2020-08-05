context("C++")
test_that("Catch unit tests pass", {
    setup_gmock()
    expect_cpp_tests_pass("malariasimulation")
})
