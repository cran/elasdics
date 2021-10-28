#' @title Retransform srv curve back to curve
#' @param srv_curve srv curve as a function of one parameter,
#' needs to be vectorised.
#' @param t time points at which the resulting curve shall be evaluated.
#' @return a \code{matrix} with curve evaluations at time points t in its columns,
#' rows correspond to coordinate directions


srvf_to_curve <- function(t, srv_curve) {
  ndim <- length(srv_curve(0))
  integrand <- function(t) sweep( srv_curve(t), 2, sqrt(colSums(srv_curve(t)^2)), "*" )
  piece_integrals <- lapply(1:ndim, function(j){
    integrand_j <- function(t) integrand(t)[j,]
    piece_integrate <- function(i, f) integrate(f = f, lower = t[i], upper = t[i+1],
                                                stop.on.error = FALSE)$value
    piece_integral <- c(0, sapply(1:(length(t)-1), piece_integrate, f = integrand_j))
  })
  piece_integrals <- do.call("rbind", piece_integrals)
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
    matrix(colMeans(data_curve[,coord_idx, drop = FALSE]), nrow = nrow(data_curve),
           ncol = ncol(data_curve[,coord_idx, drop = FALSE]), byrow = TRUE)
  data_curve
}
