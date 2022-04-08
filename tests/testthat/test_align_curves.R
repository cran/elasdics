test_that("compute_distance is zero",{
  data_curve <- data.frame(x1 = sin(1:6/3*pi), x2 = cos(1:6/3*pi))
  srv_data <- get_srv_from_points(data_curve)
  expect_equal(compute_distance(srv_data, srv_data, 0:5/5, closed = FALSE), 0,
               tolerance=1e-4)
  expect_equal(compute_distance(srv_data, srv_data, 0:5/5, closed = TRUE), 0,
               tolerance=1e-4)
})


test_that("three points shape",{
  data_curve1 <- data.frame(x1 = sin(1:6/3*pi), x2 = cos(1:6/3*pi))
  data_curve2 <- data.frame(x1 = sin(0:2*pi), x2 = cos(0:2*pi))
  expect_equal(align_curves(data_curve1, data_curve2)$data_curve2_aligned$t_optim[[2]], 0.4)
})

test_that("input checking parametrisation", {
  data_curve <- data.frame(x1 = sin(1:6), x2 = cos(1:6), t = 0:5)
  expect_error(align_curves(data_curve, data_curve),
               "Parametrisation t needs to be within 0 and 1 and increasing!")
  data_curve$t <- 1:6/6
  expect_error(align_curves(data_curve, data_curve),
               "Parametrisation t needs to start at 0!")
  data_curve$t <- 0:5/6
  expect_error(align_curves(data_curve, data_curve),
               "Last value of parametrisation t needs to be 1!")
  data_curve <- rbind(data_curve[1,], data_curve)
  data_curve$t <- NULL
  expect_warning(align_curves(data_curve, data_curve), "Duplicated points in data curves have been removed!")
})

test_that("input checking closed curves", {
  data_curve <- data.frame(x1 = sin(1:6), x2 = cos(1:6), t = 0:5/5)
  expect_error(align_curves(data_curve, data_curve, closed = TRUE),
               "Curve is not closed")
  data_curve$t <- NULL
  data_curve[nrow(data_curve),] <- data_curve[1,] + .Machine$double.eps
  expect_equal(align_curves(data_curve, data_curve, closed = TRUE)$elastic_dist, 0,
               tolerance = 10^-4)
})

test_that("it doesn't matter which coloum t is", {
  data_curve1 <- data.frame(x1 = sin(1:6), x2 = cos(1:6), t = 0:5/5)
  data_curve2 <- data.frame(t = 0:5/5, x1 = sin(1:6), x2 = cos(1:6))
  expect_equal(align_curves(data_curve1, data_curve2)$elastic_dist,0)
})

test_that("curves get closed", {
  data_curve1 <- data.frame(x1 = sin(1:5), x2 = cos(1:5))
  data_curve2 <- rbind(data_curve1, data_curve1[1,])
  data_curve2$t <- 0:5/5
  dist <- align_curves(data_curve1, data_curve2, closed = TRUE)$elastic_dist
  expect_equal(dist, 0, tolerance = 10^-4)
  data_curve1$t <- 0:4/5
  dist <- align_curves(data_curve1, data_curve2, closed = TRUE)$elastic_dist
  expect_equal(dist, 0, tolerance = 10^-4)
})

test_that("more complicated open curves", {
  data_curve1 <- data.frame(x1 = 1:20*sin(1:20), x2 = 1:10*cos(1:20))
  data_curve2 <- data.frame(x1 = sin(1:20), x2 = cos(1:20/2))
  expect_equal(align_curves(data_curve1, data_curve2)$elastic_dist, 10.062,
               tolerance = 1e-2)
})

test_that("3d curves", {
  data_curve1 <- data.frame(x1 = 1:20*sin(1:20), x2 = 1:10*cos(1:20), x3 = tan(1:20))
  data_curve2 <- data.frame(x1 = sin(1:20), x2 = cos(1:20/2), x3 = atan(1:20))
  expect_equal(align_curves(data_curve1, data_curve2)$elastic_dist, 24.38,
               tolerance = 1e-2)
})

test_that("1d curves", {
  data_curve1 <- data.frame(x1 = 1:20*sin(1:20)/20)
  data_curve2 <- data.frame(x1 = cos(1:15))
  expect_warning(align_curves(data_curve1, data_curve2))
  data_curve1$t <- 0:(nrow(data_curve1) - 1)/(nrow(data_curve1) - 1)
  expect_warning(align_curves(data_curve1, data_curve2))
})

test_that("same dim for both curves", {
  data_curve1 <- data.frame(x1 = 1:20*sin(1:20), x2 = 1:10*cos(1:20), x3 = tan(1:20))
  data_curve2 <- data.frame(x1 = sin(1:20), x2 = cos(1:20/2))
  expect_error(align_curves(data_curve1, data_curve2), "Both curves must have same number of dimensions!")
})

test_that("closed curve with points mapped to same optimal time", {
  t_grid <- seq(0,1, length = 10)
  data_curve1 <- data.frame(x1 = t_grid, x2 = sin(9*t_grid))
  data_curve2 <- data.frame(x1 = t_grid, x2 = sin(3*t_grid))
  expect_equal(sum(diff(align_curves(data_curve2, data_curve1, closed = TRUE)$data_curve2_aligned$t_optim) == 0) > 0, TRUE)
})

