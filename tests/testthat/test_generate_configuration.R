context("Check generate_configuration")

test_that("Test if default output is in desired format", {

  #check default configuration
  default_configuration <- generate_configuration()
  expect_is(default_configuration, "list")
  expect_is(default_configuration$n_iter, "numeric")
  expect_is(default_configuration$n_burnin, "numeric")
  expect_is(default_configuration$set_seed, "logical")

})


