#' @title Retransform srv curve back to curve
#' @param srv_curve srv curve as a function of one parameter,
#' needs to be vectorised.
#' @param t time points at which the resulting curve shall be evaluated.
#' @return a \code{matrix} with curve evaluations at time points t in its columns,
#' rows correspond to coordinate directions


srvf_to_curve <- function(t, srv_curve) {
  integrand <- function(t) sweep( srv_curve(t), 2, sqrt(colSums(srv_curve(t)^2)), "*" )
  integrand_1 <- function(t) integrand(t)[1,]
  integrand_2 <- function(t) integrand(t)[2,]
  piece_integrate <- function(i, f) integrate(f = f, lower = t[i], upper = t[i+1],
                                              stop.on.error = FALSE)$value
  piece_integrals <- rbind( sapply(1:(length(t)-1), piece_integrate, f = integrand_1),
                            sapply(1:(length(t)-1), piece_integrate, f = integrand_2))
  piece_integrals <- cbind(c(0,0), piece_integrals)
  t(apply(piece_integrals, 1, cumsum))
}


#' @title Centers curves for plotting
#' @param data_curve curve data
#' @return a \code{data.frame} with evaluations of the curve
#' centered at the origin
#' @export

center_curve <- function(data_curve){
  coord_idx <- !(colnames(data_curve) %in% c("t", "id"))
  data_curve[,coord_idx] <- data_curve[,coord_idx] -
    matrix(colMeans(data_curve[,coord_idx]), nrow = nrow(data_curve),
           ncol = ncol(data_curve[,coord_idx]), byrow = TRUE)
  data_curve
}
