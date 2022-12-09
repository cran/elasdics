test_that("input checkin",{
  data_curve1 <- data.frame("y1" = sin(1:7/4*pi), "y2" = cos(1:7/4*pi), "y3" = 1)
  data_curve2 <- data.frame("y1" = sin(1:15/8*pi), "y2" = cos(1:15/8*pi), "y3" = 1)
  data_curves <- list(data_curve1, data_curve2)
  reg_model <- fit_elastic_regression(data_curves ~ x, data_curves = data_curves,
                                      x_data = data.frame("x" = c(-1,1)))
  expect_error(plot(reg_model), "Plotting option only for planar curves!")
})

test_that("correct number of predicted curves",{
  data_curve1 <- data.frame("y1" = sin(1:7/4*pi), "y2" = cos(1:7/4*pi))
  data_curve2 <- data.frame("y1" = sin(1:15/8*pi), "y2" = cos(1:15/8*pi))
  data_curves <- list(data_curve1, data_curve2)
  reg_model <- fit_elastic_regression(data_curves ~ x, data_curves = data_curves,
                                      x_data = data.frame("x" = c(-1,1)))
  expect_equal(length(predict(reg_model)), nrow(reg_model$x_data))
  expect_warning(plot(reg_model), regexp = NA)
})

test_that("predict closed polygon model",{
  data_curve1 <- data.frame(x1 = sin(1:6/4*pi), x2 = cos(1:6/4*pi))
  data_curve2 <- data_curve <- data.frame(x1 = sin(1:15/8*pi), x2 = cos(1:15/8*pi))
  data_curves <- list(data_curve1, data_curve2)
  reg_model <- fit_elastic_regression(data_curves ~ x, data_curves = data_curves, type = "polygon",
                                      knots = seq(0,1,0.02),
                                      x_data = data.frame("x" = c(-1,1)), closed =TRUE)
  expect_equal(length(predict(reg_model)), nrow(reg_model$x_data))
})

test_that("intercept model only",{
  data_curve1 <- data.frame("y1" = sin(1:7/4*pi), "y2" = cos(1:7/4*pi))
  data_curve2 <- data.frame("y1" = sin(1:15/8*pi), "y2" = cos(1:15/8*pi))
  data_curves <- list(data_curve1, data_curve2)
  reg_model <- fit_elastic_regression(data_curves ~ 1, data_curves = data_curves,
                                      x_data = data.frame("x" = c(-1,1)))
  expect_equal(length(predict(reg_model)), 1)
  expect_warning(plot(reg_model), regexp = NA)
})

test_that("plot factor variables",{
  data_curve1 <- data.frame("y1" = sin(1:7/4*pi), "y2" = cos(1:7/4*pi))
  data_curve2 <- data.frame("y1" = sin(1:15/8*pi), "y2" = cos(1:15/8*pi))
  data_curves <- list(data_curve1, data_curve2)
  reg_model <- fit_elastic_regression(data_curves ~ x, data_curves = data_curves,
                                      x_data = data.frame("x" = factor(c(-1,1))))
  expect_warning(plot(reg_model), regexp = NA)
})
