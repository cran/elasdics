#' @title Close open curve via projection on derivative level.
#' @param data_curve \code{data.frame} with values of the curve.
#' @return a \code{data.frame} with closed curve.

project_curve_on_closed <- function(data_curve){
  delta <- data_curve[nrow(data_curve), ] - data_curve[1, ]
  delta_parts <- sapply(delta, seq, from = 0, length = nrow(data_curve))
  data_curve <- data.frame(data_curve - delta_parts)

  #center curve
  offset <- apply(get_evals(cbind("t" = get_arc_length_param(data_curve),
                                  data_curve)), 2, mean)
  data.frame(t(t(data_curve) - offset))
}
