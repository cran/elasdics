test_that("input checking",{
  data_curve1 <- data.frame(x1 = sin(1:7/4*pi), x2 = cos(1:7/4*pi))
  data_curve2 <- data_curve <- data.frame(x1 = sin(1:15/8*pi), x2 = cos(1:15/8*pi))
  data_curves <- list(data_curve1, data_curve2)
  expect_error(fit_elastic_regression(y ~ x, data_curves = data_curves,
                                      x_data = data.frame("x" = c(-1,1))),
               "formula must be of form data_curves ~ ...")
})

test_that("correct number of estimated coefficients",{
  data_curve1 <- data.frame(x1 = sin(1:7/4*pi), x2 = cos(1:7/4*pi))
  data_curve2 <- data_curve <- data.frame(x1 = sin(1:15/8*pi), x2 = cos(1:15/8*pi))
  data_curves <- list(data_curve1, data_curve2)
  reg_model <- fit_elastic_regression(data_curves ~ x, data_curves = data_curves, x_data = data.frame("x" = c(-1,1)))
  expect_equal(names(reg_model$coefs_list), c("beta_0", "beta_1"))
})

test_that("closed curves",{
  data_curve1 <- data.frame(x1 = sin(1:6/4*pi), x2 = cos(1:6/4*pi))
  data_curve2 <- data_curve <- data.frame(x1 = sin(1:15/8*pi), x2 = cos(1:15/8*pi))
  data_curves <- list(data_curve1, data_curve2)
  reg_model <- fit_elastic_regression(data_curves ~ x, data_curves = data_curves, type = "polygon",
                                      knots = seq(0,1,0.1),
                                      x_data = data.frame("x" = c(-1,1)), closed =TRUE)
  expect_equal(names(reg_model$coefs_list), c("beta_0", "beta_1"))
  reg_model <- fit_elastic_regression(data_curves ~ x, data_curves = data_curves,
                                      x_data = data.frame("x" = c(-1,1)), closed =TRUE,
                                      max_iter = 1)
  expect_equal(names(reg_model$coefs_list), c("beta_0", "beta_1"))
})
