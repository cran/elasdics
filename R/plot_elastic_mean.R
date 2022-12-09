#' Plot method for planar elastic mean curves
#' @description Plots objects of class \code{elastic_mean}.
#' @param x object of class \code{elastic_mean},
#' usually a result of a call to \code{\link{compute_elastic_mean}}
#' @param asp numeric, giving the aspect ratio of the two coordinates,
#' see \code{\link{plot.window}} for details.
#' @param col color of the mean curve.
#' @param ... further plotting parameters.
#' @importFrom graphics plot lines
#' @return No value
#' @export
#'
#' @seealso For examples see documentation of \code{\link{compute_elastic_mean}}.

plot.elastic_mean <- function(x, asp = 1, col = "red", ...){
  if(!(ncol(x$coefs) %in% 1:2)){
    stop("Plotting option only for functions and planar curves!")
  }
  if(ncol(x$coefs) == 1){
    data_curves <- lapply(x$data_curves, function(data) data[, -1])
    data_range <- range(do.call("c", data_curves))
    #empty plot
    plot(NULL, xlim = c(0,1), ylim = data_range, xlab = "t",
         ylab = colnames(data_curves[[1]])[2])
    #plot data
    lapply(data_curves, lines, col = "gray")
    #plot mean
    t_grid <- seq(0,1, 0.01)
    mean_data <- data.frame("t" = t_grid, get_evals(x, t_grid = t_grid))
    lines(mean_data, col = col, lwd = 2)

  } else {
    data_curves <- lapply(x$data_curves, center_curve)
    data_curves <- lapply(data_curves, function(data) data[,colnames(x$coefs)])
    data_all <- do.call("rbind", data_curves)

    #empty plot
    plot(NULL, xlim = range(data_all[,1]), ylim = range(data_all[,2]), xlab = colnames(x$coefs)[1],
         ylab = colnames(x$coefs)[2], asp = asp)
    #plot data
    invisible(lapply(data_curves, lines, col = "gray"))

    #plot mean
    lines(get_evals(x), col = col, lwd = 2)
  }
}
