test_that("open curve",{
  data_curve1 <- data.frame(x1 = c(1, 0.5, -1, -1), x2 = c(1, -0.5, -1, 1))
  data_curve2 <- data.frame(x1 = c(0.1,0.7)*sin(1:6), x2 = cos(1:6))
  aligned_curves <- align_curves(data_curve1, data_curve2)
  expect_warning(plot(aligned_curves), regexp = NA)
})

test_that("closed curves",{
  data_curve1 <- data.frame(x1 = sin(0:12/5), x2 = cos(0:12/5))
  data_curve2 <- data.frame(x1 = c(1, 0.5, -1, -1), x2 = c(1, -0.5, -1, 1))
  aligned_curves_closed <- align_curves(data_curve1, data_curve2, closed = TRUE)
  expect_warning(plot(aligned_curves_closed), regexp = NA)

  data_curve1 <- data.frame(x1 = sin(0:4/5), x2 = cos(0:4/5))
  data_curve2 <- data.frame(x1 = c(1, 0.5, -1, -1), x2 = c(1, -0.5, -1, 1))
  aligned_curves_closed <- align_curves(data_curve2, data_curve1, closed = TRUE)
  expect_warning(plot(aligned_curves_closed), regexp = NA)
})

test_that("1d function",{
  t_grid = 0:20/20
  data_curve1 <- data.frame(t = t_grid, x = sin(6*t_grid))
  data_curve2 <- data.frame(t = t_grid, x = cos(3*t_grid - 1))
  expect_warning(aligned_curves <- align_curves(data_curve1, data_curve2))
  expect_warning(plot(aligned_curves), regexp = NA)
})

test_that("get_warping",{
  data_curve1 <- data.frame(x1 = c(1, 0.5, -1, -1), x2 = c(1, -0.5, -1, 1))
  data_curve2 <- data.frame(x1 = c(0.1,0.7)*sin(3:8), x2 = cos(3:8))
  aligned_curves <- align_curves(data_curve1, data_curve2)
  expect_true(all(aligned_curves$data_curve2_aligned$t_optim %in% get_warping(aligned_curves)$t))
})

test_that("only for 1d or 2d curves",{
  data_curve1 <- data.frame(x1 = c(1, 0.5, -1, -1), x2 = c(1, -0.5, -1, 1), x3 = c(1, -0.5, -1, 1))
  data_curve2 <- data.frame(x1 = c(0.1,0.7)*sin(1:6), x2 = cos(1:6), x3 = cos(1:6))
  aligned_curves <- align_curves(data_curve1, data_curve2)
  expect_error(plot(aligned_curves), "Plotting method only for 1d functions and 2d curves available!")
})
