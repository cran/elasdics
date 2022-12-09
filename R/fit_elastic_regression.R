#' Compute a elastic mean for a collection of curves
#' @name fit_elastic_regression
#' @description Computes a Fréchet mean for the curves stored in \code{data_curves} with respect
#' to the elastic distance. Constructor function for class \code{elastic_reg_model}.
#' @param formula an object of class "formula" of the form data_curves ~ ...".
#' @param data_curves list of \code{data.frame}s with observed points in each row. Each
#' variable is one coordinate direction. If there is a variable \code{t},
#' it is treated as the time parametrization, not as an additional coordinate.
#' @param x_data a \code{data.frame} with covariates.
#' @param knots set of knots for the parameter curves of the regression model
#' @param type if "smooth" linear srv-splines are used which results in a differentiable mean curve
#' if "polygon" the mean will be piecewise linear.
#' @param closed \code{TRUE} if the curves should be treated as closed.
#' @param eps the algorithm stops if L2 norm of coefficients changes less
#' @param max_iter maximal number of iterations
#' @return an object of class \code{elastic_reg_model}, which is a \code{list}
#' with entries
#'   \item{type}{"smooth" if linear srv-splines or
#'   "polygon" if constant srv-splines were used}
#'   \item{coefs}{spline coeffiecients}
#'   \item{knots}{spline knots}
#'   \item{data_curves}{list of \code{data.frame}s with observed points in each row.
#'   First variable \code{t} gives the initial parametrization, second variable \code{t_optim}
#'   the optimal parametrization when the curve is aligned to the model prediction.}
#'   \item{closed}{\code{TRUE} if the regression model fitted closed curves.}
#' @export
#' @exportClass elastic_reg_model
#' @importFrom splines splineDesign
#' @examples
#' curve <- function(x_1, x_2, t){
#'   rbind(2*t*cos(6*t) - x_1*t , x_2*t*sin(6*t))
#' }
#' set.seed(18)
#' x_data <- data.frame(x_1 = runif(10,-1,1), x_2 = runif(10,-1,1))
#' data_curves <- apply(x_data, 1, function(x){
#'   m <- sample(10:15, 1)
#'   delta <- abs(rnorm(m, mean = 1, sd = 0.05))
#'   t <- cumsum(delta)/sum(delta)
#'   data.frame(t(curve((x[1] + 1), (x[2] + 2), t))
#'    + 0.07*t*matrix(cumsum(rnorm(2*length(delta))), ncol = 2))
#' })
#' reg_model <- fit_elastic_regression(data_curves ~ x_1 + x_2,
#'                                     data_curves = data_curves, x_data = x_data)
#' plot(reg_model)

fit_elastic_regression <- function(formula, data_curves, x_data, knots = seq(0,1,0.2), type = "smooth",
                                   closed = FALSE, max_iter = 10, eps = 0.001){
  formula <- as.formula(formula)
  if(formula[[2]] != "data_curves") stop("formula must be of form data_curves ~ ...")

  # create x dataframe
  x_data <- data.frame(x_data)

  # create x model matrix
  x_model_matrix <- model.matrix(formula, cbind("data_curves" = 1, x_data))

  # input checking for closed curves
  if(closed) data_curves <- lapply(data_curves, check_closed)

  #remove duplicated points
  data_curves <- lapply(data_curves, remove_duplicate, closed = closed)

  #compute srv
  srv_data <- lapply(data_curves, get_srv_from_points)

  #initial alignment as optimal alignment to the mean
  mean <- compute_elastic_mean(data_curves, knots = seq(0,1,0.01), max_iter = 100, type = "polygon",
                               eps = 0.001, closed = closed)

  data_curves <- lapply(mean$data_curves, function(data_curve){
    if(data_curve$t_optim[nrow(data_curve)] == 0){
      data_curve$t_optim[nrow(data_curve)] <- 1
    }
    attr(data_curve, "dist_to_prediction") <- attributes(data_curve)$dist_to_mean
    attr(data_curve, "dist_to_mean") <- NULL
    data_curve
  })

  srv_data_initial <- lapply(1:length(srv_data), function(i){
    dat <- srv_data[[i]]
    idx_start <- which(data_curves[[i]]$t_optim == 0)[1]:nrow(dat)
    dat <- rbind(dat[idx_start,], dat[-idx_start,])
    dat$t <- dat$t - dat$t[1]
    dat$t <- sapply(dat$t, function(t) ifelse(t < 0, t + 1, t))
    dat
  })
  data_curves_t <- lapply(srv_data_initial, function(srv_data_initial){
    cbind("t" = c(srv_data_initial$t, 1), get_points_from_srv(srv_data_initial))
  })

  #initiate values
  t_optims <- lapply(mean$data_curves, function(curve){
    c(sort(curve$t_optim[-length(curve$t_optim)]), 1)
  })
  coefs_list <- 0

  for(i in 1:max_iter){
    coefs_list_old <- coefs_list

    if(type == "polygon"){
      srv_data_discrete <- lapply(1:length(data_curves), function(i){
        data_curves_t[[i]]$t <- t_optims[[i]]
        curve_data <- get_evals(data_curves_t[[i]], t_grid = knots)
        unlist(get_srv_from_points(cbind("t" = knots, curve_data))[, -1])
      })
      coefs_discrete <- sapply(1:length(srv_data_discrete[[1]]), function(i){
        response <- sapply(srv_data_discrete, '[[', i)
        lm(response ~ -1 + x_model_matrix)$coefficients
      })
      coefs_list <- lapply(1:nrow(coefs_discrete), function(i){
        coefs <- matrix(coefs_discrete[i,], ncol = ncol(data_curves[[1]]) -2)
        colnames(coefs) <- colnames(srv_data[[1]][,-1])
        coefs
      })

    } else {
      data_curves_t_now <- lapply(1:length(data_curves_t), function(i){
        idx <- c(TRUE, diff(t_optims[[i]]) != 0)
        data.frame("t" = t_optims[[i]][idx], data_curves_t[[i]][idx, -1])
      })

      model_data <- get_reg_model_data(data_curves_t_now, seq(0,1,0.01), x_model_matrix)
      design_mat <- make_reg_design(model_data[, 1:ncol(x_model_matrix)], model_data$m_long,
                                    knots = knots, type = type, closed = closed)

      coefs_long <- apply(model_data[,-(1:(ncol(x_model_matrix) + 1)), drop = FALSE], 2, function(q_m_x_long){
        q_m_x_long[!is.finite(q_m_x_long)] <- NA
        coef(lm(q_m_x_long ~ -1 + design_mat))
      })

      n_coefs <- nrow(coefs_long)/ncol(x_model_matrix)
      coefs_list <- lapply(1:ncol(x_model_matrix), function(i){
        coefs <- coefs_long[n_coefs*(i - 1) + 1:n_coefs, ]
        colnames(coefs) <- colnames(srv_data[[1]][,-1])
        rownames(coefs) <- paste0("coef_", 1:n_coefs)
        coefs
      })
    }
    names(coefs_list) <- paste0("beta_", 1:ncol(x_model_matrix) - 1)

    #stop if coefficients don't change much anymore
    stop_crit <- sum((unlist(coefs_list) - unlist(coefs_list_old))^2)/sum(unlist(coefs_list)^2)
    if(stop_crit < eps | max_iter == 0){
      elastic_reg_model <- list("coefs_list" = coefs_list, "formula" = formula,
                                "data_curves" = data_curves, "x_data" = x_data,
                                "knots" = knots, "type" = type, "closed" = closed)
      class(elastic_reg_model) <- "elastic_reg_model"
      return(elastic_reg_model)
    }
    ############################################################################
    # update warping alignment
    t_optims <- lapply(1:length(data_curves), function(i){
      x_i <- x_model_matrix[i,]
      pred_coefs <- Reduce("+", lapply(1:length(x_i), function(j) x_i[j]*coefs_list[[j]]))

      if(type == "smooth"){
        pfun <- function(t){
          t(make_reg_design(x_design_mat = NULL, t = t, knots = knots,
                            type = type, closed = closed) %*% pred_coefs)
        }

        find_optimal_t(srv_curve = pfun, s = c(srv_data_initial[[i]][,1],1),
                       q = t(srv_data_initial[[i]][,-1]), initial_t = t_optims[[i]],
                       eps = 0.01)
      } else {
        find_optimal_t_discrete(r = knots, p = t(pred_coefs), s = c(srv_data_initial[[i]][,1],1),
                                q = t(srv_data_initial[[i]][,-1]), initial_t = t_optims[[i]],
                                eps = 0.01)
      }
    })
    # add optimal time parametrisation to data curves
    for(i in 1:length(data_curves)){
      warping <- data_curves[[i]][-nrow(data_curves[[i]]),1:2]
      warping <- warping[order(warping$t_optim), ]
      warping$t_optim <- t_optims[[i]][-length(t_optims[[i]])]
      warping <- warping[order(warping$t), ]
      data_curves[[i]][,1:2] <- rbind(warping, c(1, ifelse(warping[1,2] == 0, 1, warping[1,2])))
      attr(data_curves[[i]], "dist_to_prediction") <- attributes(t_optims[[i]])$dist
    }
  }

  elastic_reg_model <- list("coefs_list" = coefs_list, "formula" = formula,
                            "data_curves" = data_curves, "x_data" = x_data,
                            "knots" = knots, "type" = type, "closed" = closed)
  class(elastic_reg_model) <- "elastic_reg_model"
  return(elastic_reg_model)
}

get_reg_model_data <- function(data_curves_t, t_grid, x_model_data){
  q_m <- lapply(data_curves_t, function(data_curve){
    data_curve <- data.frame("t" = t_grid, get_evals(data_curve, t_grid = t_grid))
    get_srv_from_points(data_curve)[,-1]
  })

  #convert in long format
  x_long <- do.call("rbind", lapply(1:nrow(x_model_data), function(i){
    sapply(x_model_data[i,], function(x) rep(x, length(t_grid) - 1))
    }))
  m_long <- rep(t_grid[-length(t_grid)] + diff(t_grid)/2, length(q_m))
  q_m_long <- do.call(rbind, q_m)
  data.frame("x_long" = x_long, "m_long" = m_long, "q_m_long" = q_m_long)
}

# creating the design matrix
make_reg_design <- function(x_design_mat, t, knots, type = "smooth", closed = FALSE) {
  deg <- ifelse(type == "smooth", 1, 0)
  spline_design_mat <- splines::splineDesign(knots = c(rep(0, deg), knots, rep(1, deg)),
                                    x = t, outer.ok = TRUE, ord = deg + 1)
  if(closed == TRUE & type == "smooth"){
    spline_design_mat[,1] <- spline_design_mat[,1] + spline_design_mat[, ncol(spline_design_mat)]
    spline_design_mat <- spline_design_mat[,-ncol(spline_design_mat)]
  }

  if(is.null(x_design_mat)) return(spline_design_mat)
  x_design_mat <- as.matrix(x_design_mat)
  t(sapply(1:nrow(x_design_mat), function(i) kronecker(x_design_mat[i,], spline_design_mat[i,])))
}

