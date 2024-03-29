#' Align two curves measured at discrete points
#' @name align_curves
#' @description Finds the optimal reparametrization of the second curve (stored in
#' \code{data_curve2}) to the first one (stored in \code{data_curve1}) with respect
#' to the elastic distance. Constructor function for class \code{aligned_curves}.
#' @param data_curve1 \code{data.frame} with observed points in each row. Each
#' variable is one coordinate direction. If there is a variable \code{t},
#' it is treated as the time parametrization, not as an additional coordinate.
#' @param data_curve2 same as \code{data_curve1}
#' @param closed \code{TRUE} if the curves should be treated as closed.
#' @param eps convergence tolerance
#' @return an object of class \code{aligned_curves}, which is a \code{list}
#' with entries
#'   \item{data_curve1}{\code{data_curve1} with parametrization variable \code{t}}
#'   \item{data_curve2_aligned}{\code{data_curve2} with initial parametrization
#'   variable \code{t} and optimal parametrization \code{t_optim}}
#'   \item{elastic_dist}{elastic distance between curve1 and curve2}
#'   \item{closed}{\code{TRUE} if the curves should have been treated as closed.}
#' @export
#' @examples
#' #open curves
#' data_curve1 <- data.frame(x1 = c(1, 0.5, -1, -1), x2 = c(1, -0.5, -1, 1))
#' data_curve2 <- data.frame(x1 = c(0.1,0.7)*sin(1:6), x2 = cos(1:6))
#' aligned_curves <- align_curves(data_curve1, data_curve2)
#' plot(aligned_curves)
#'
#' #different parametrization of the first curve
#' data_curve1$t <- 0:3/3
#' align_curves(data_curve1, data_curve2)
#'
#' #closed curves
#' data_curve1 <- data.frame(x1 = sin(0:12/5), x2 = cos(0:12/5))
#' data_curve2 <- data.frame(x1 = c(1, 0.5, -1, -1), x2 = c(1, -0.5, -1, 1))
#' aligned_curves_closed <- align_curves(data_curve1, data_curve2, closed = TRUE)
#' plot(aligned_curves_closed, asp = 1)


align_curves <- function(data_curve1, data_curve2, closed = FALSE, eps = 0.01){
  # input checking given parametrization t
  if("t" %in% names(data_curve1)) check_param(data_curve1, closed)
  if("t" %in% names(data_curve2)) check_param(data_curve2, closed)

  #remove duplicated points
  data_curve1 <- remove_duplicate(data_curve1, closed = closed)
  data_curve2 <- remove_duplicate(data_curve2, closed = closed)
  if(attributes(data_curve1)$points_rm | attributes(data_curve2)$points_rm){
    warning("Duplicated points in data curves have been removed!")
  }

  # input checking for closed curves
  if(closed){
    data_curve1 <- check_closed(data_curve1)
    data_curve2 <- check_closed(data_curve2)
  }

  srv_data_1 <- get_srv_from_points(data_curve1)
  srv_data_2 <- get_srv_from_points(data_curve2)
  if(ncol(srv_data_1) != ncol(srv_data_2)) stop("Both curves must have same number of dimensions!")
  #remove parametrization
  if(ncol(srv_data_1) == 2) warning("This package was designed to analyse curve data in d-dimensions with d > 1.
                                    It might still be used for functional data (d = 1) but results might not be satisfing")
  if("t" %in% names(data_curve1)) data_curve1$t  <- NULL
  if("t" %in% names(data_curve2)) data_curve2$t  <- NULL
  # after computing the srv transformations the parametrization t is definitly
  # in the first column
  if(closed){
    #pre alignment
    t <- get_arc_length_param(data_curve2)
    loss <- sapply(1:length(srv_data_2$t), function(i){
      t_new <- t - t[i]
      compute_distance(srv_data_1, srv_data_2, t_new, closed)
    })
    initial_t <- t - t[which.min(loss)]
    optimal_t <- find_optimal_t_discrete_closed(c(srv_data_1$t,1), t(srv_data_1[,-1]),
                                         c(srv_data_2$t,1), t(srv_data_2[,-1]),
                                         initial_t = initial_t, eps = eps)
    t_optim <- optimal_t + (optimal_t < 0)
    t_optim <- ifelse(t_optim > 1, t_optim - 1, t_optim)
  } else {
    initial_t <- get_arc_length_param(data_curve2)
    t_optim <- find_optimal_t_discrete(c(srv_data_1$t,1), t(srv_data_1[,-1]),
                                       c(srv_data_2$t,1), t(srv_data_2[,-1]),
                                       initial_t = initial_t, eps = eps)
  }

  elastic_dist <- compute_distance(srv_data_1, srv_data_2, t_optim, closed)

  data_curve1 <- cbind(t = c(srv_data_1$t, 1), data_curve1)
  data_curve2_aligned <- cbind(t = c(srv_data_2$t, 1),
                               t_optim = t_optim, data_curve2)

  aligned_curves <- list("data_curve1" = data_curve1,
                         "data_curve2_aligned" = data_curve2_aligned,
                         "elastic_dist" = elastic_dist,
                         "closed" = closed)
  class(aligned_curves) <- "aligned_curves"
  return(aligned_curves)

}

#' Input checking for given parametrization
#' @inheritParams align_curves
#' @param data_curve data of curve like in \code{align_curves}
#' @noRd

check_param <- function(data_curve, closed = closed){
  if(!(data_curve$t[1] >= 0 & data_curve$t[nrow(data_curve)] <= 1 &
       all(diff(data_curve$t) >= 0))){
    stop("Parametrization t needs to be within 0 and 1 and increasing!")
    }
  if(data_curve$t[1] != 0){
    stop("Parametrization t needs to start at 0!")
    }
  if(!closed & data_curve$t[nrow(data_curve)] != 1){
    stop("Last value of parametrization t needs to be 1!")
  }
}



#' Input checking for closed curves
#' @param data_curve data of curve like in \code{align_curves}
#' @noRd

check_closed <- function(data_curve){
  if("t" %in% names(data_curve)){
    if(data_curve$t[nrow(data_curve)] == 1){
      if(!all(data_curve[1, names(data_curve) != "t"] ==
         data_curve[nrow(data_curve), names(data_curve) != "t"])){
        stop("Curve is not closed")
      }
    } else {
      data_curve <- rbind(data_curve, data_curve[1,])
      data_curve$t[nrow(data_curve)] <- 1
    }
  } else {
    err_non_closed <- sum((data_curve[1,] - data_curve[nrow(data_curve),])^2)/sum(data_curve[1,]^2)
    if(err_non_closed > sqrt(.Machine$double.eps)) {
      data_curve <- rbind(data_curve, data_curve[1,])
    } else {
      data_curve[nrow(data_curve),] <- data_curve[1,]
    }
  }
  return(data_curve)
}

#' Remove duplicated points in data curves
#' @param data_curve data of curve like in \code{align_curves}
#' @noRd

remove_duplicate <- function(data_curve, closed){
  if(ncol(data_curve) == 1){
    attr(data_curve, "points_rm") <- FALSE
    return(data_curve)
  }
  points <- as.data.frame(data_curve)
  try(points$t <- NULL, silent = TRUE)
  moves <- c(TRUE, rowSums(apply(points, 2, diff)^2) > max(points)*.Machine$double.eps)
  data_curve <- data_curve[moves,]
  attr(data_curve, "points_rm") <- !all(moves)
  if(!is.null(data_curve$t) & !closed) data_curve$t[nrow(data_curve)] <- 1
  data_curve
}

#' Computes elastic distance
#' @noRd
#' @param srv_data_1 srv of curve1
#' @param srv_data_2 srv of curve2
#' @param t_optim optimal parametrization of curve2

compute_distance<- function(srv_data_1, srv_data_2, t_optim, closed){
  norm_1 <- sum(t(srv_data_1[,-1]^2)%*%diff(c(srv_data_1$t,1)))
  norm_2 <- sum(t(srv_data_2[,-1]^2)%*%diff(c(srv_data_2$t,1)))
  if(closed){
    srv_data_1_extended <- rbind(srv_data_1, srv_data_1, srv_data_1)
    srv_data_1_extended$t <- c(srv_data_1$t - 1, srv_data_1$t, srv_data_1$t + 1)
    cross_prod <- get_loss_discrete(t = t_optim, srv_data_1 = srv_data_1_extended, srv_data_2)
  } else {
    cross_prod <- get_loss_discrete(t = t_optim, srv_data_1, srv_data_2)
  }
  dist_squared <- norm_1 + norm_2 - 2*cross_prod
  sqrt(ifelse(dist_squared <= 0, 0, dist_squared))
}


