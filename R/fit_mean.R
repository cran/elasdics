#' Fitting function for open curves
#' @name fit_mean
#' @description Fits an elastic mean for open curves. Is usually called from
#' \code{\link{compute_elastic_mean}}.
#' @param srv_data_curves list of \code{data.frame}s with srv vectors in each row.
#' Usually a result of a call to \code{\link{get_srv_from_points}}
#' @param knots set of knots for the mean spline curve
#' @param max_iter maximal number of iterations
#' @param eps the algorithm stops if L2 norm of coefficients changes less
#' @param type if "smooth" linear srv-splines are used which results in a differentiable mean curve
#' if "polygon" the mean will be piecewise linear.
#' @return a \code{list}
#' with entries
#'   \item{type}{"smooth" or "polygon"}
#'   \item{coefs}{\code{coefs} srv spline coefficients of the estimated mean}
#'   \item{knots}{spline knots}
#'   \item{t_optims}{optimal parametrisation}

fit_mean <- function(srv_data_curves, knots, max_iter, type, eps){
  #initial param and coefs
  t_optims <- lapply(srv_data_curves, function(srv_data_curve){
    c(srv_data_curve$t, 1)
  })
  coefs <- 0

  #iterate L_2 mean fit and warping
  for (i in 1:max_iter){
    model_data <- get_model_data(t_optims, srv_data_curves, knots, type)

    # fit model
    coefs_old <- coefs
    coefs <- apply(model_data[,-1, drop = FALSE], 2, function(q_m_x_long){
      q_m_x_long[!is.finite(q_m_x_long)] <- NA
      coef(lm(q_m_x_long ~ -1 + make_design(model_data[,1], knots = knots,
                                            closed = FALSE, type = type)))
    })

    #stop if coefficients don't change much anymore
    stop_crit <- sum((coefs - coefs_old)^2)/sum(coefs^2)
    if(stop_crit < eps | max_iter == 0){
      rownames(coefs) <- NULL
      colnames(coefs) <- colnames(srv_data_curves[[1]][,-1])
      return(list("type" = type, "coefs" = coefs, "knots" = knots, "t_optims" = t_optims))
    }
    ## warping fit, update t_optims
    if(type == "smooth"){
      pfun <- function(t){
        t(make_design(t, knots = knots,
                      closed = FALSE, type = type) %*% coefs)
      }

      t_optims <- lapply(1:length(srv_data_curves), function(j){
        t_optim <- find_optimal_t(srv_curve = pfun,
                                  s = c(srv_data_curves[[j]]$t, 1),
                                  q = t(srv_data_curves[[j]][,-1]),
                                  initial_t = t_optims[[j]],
                                  eps = eps*100/i)
        attr(t_optim, "dist_to_mean") <- attr(t_optim, "dist")
        attr(t_optim, "dist") <- NULL
        t_optim
      })
    } else {
      t_optims <- lapply(1:length(srv_data_curves), function(j){
        t_optim <- find_optimal_t_discrete(r = knots,
                                           p = t(coefs),
                                           s = c(srv_data_curves[[j]]$t, 1),
                                           q = t(srv_data_curves[[j]][,-1]),
                                           initial_t = t_optims[[j]],
                                           eps = eps*100/i)
        attr(t_optim, "dist_to_mean") <- attr(t_optim, "dist")
        attr(t_optim, "dist") <- NULL
        t_optim
      })
    }

  }
  warning("Stopping criteria eps has not been reached! Consider more iterations max_iter")
  rownames(coefs) <- NULL
  colnames(coefs) <- colnames(srv_data_curves[[1]][,-1])
  return(list("type" = type, "coefs" = coefs, "knots" = knots, "t_optims" = t_optims))
}

get_model_data <- function(t_optims, srv_data_curves, knots, type){
  if(type == "polygon"){
    q_m <- lapply(1: length(srv_data_curves), function(j){
      curve <- data.frame("t" = t_optims[[j]],
        get_points_from_srv(srv_data_curves[[j]]))
      curve_at_knots <- cbind("t" = knots, get_evals(curve, t_grid = knots))
      get_srv_from_points(curve_at_knots)[,-1, drop = FALSE]
    })
    m <- lapply(srv_data_curves, function(x){
      knots[-1] - 0.5*diff(knots)
    })
  } else {
    #compute warped srv vectors
    q_m <- lapply(1:length(srv_data_curves), function(j){
      old_diff <- diff(c(srv_data_curves[[j]]$t, 1))
      new_diff <- diff(t_optims[[j]])
      as.matrix(srv_data_curves[[j]][,-1])*sqrt(old_diff/new_diff)
    })

    m <- lapply(t_optims, function(t_optim){
      t_optim[-1] - 0.5*diff(t_optim)
    })
  }

  #convert in long format
  m_long <- do.call(c, m)
  q_m_long <- do.call(rbind, q_m)
  data.frame("m_long" = m_long, "q_m_long" = q_m_long)
}

# creating the design matrix
make_design <- function(t, knots, closed = FALSE, type = "smooth") {
  deg <- ifelse(type == "smooth", 1, 0)
  design_mat <- splineDesign(knots = c(rep(0, deg), knots, rep(1, deg)),
                             x = t, outer.ok = TRUE, ord = deg + 1)
  if(closed == TRUE & type == "smooth"){
    design_mat[,1] <- design_mat[,1] + design_mat[, ncol(design_mat)]
    design_mat <- design_mat[,-ncol(design_mat)]
  }
  design_mat
}
