#' Predict method for elastic regression models
#' @description predicted curves for elastic regression model objects.
#' @param object object of class \code{elastic_reg_model},
#' usually a result of a call to \code{\link{fit_elastic_regression}}
#' @param newdata an optional \code{data.frame} in which to look for variables with which to predict.
#' If not given, the fitted values are used.
#' @param t_grid grid on which the predicted curves are evaluated.
#' @param ... further arguments passed to or from other methods.
#' @return a \code{list} of \code{data.frame}s with predicted curves
#' @export
#'
#' @seealso For examples see documentation of \code{\link{fit_elastic_regression}}.

predict.elastic_reg_model <- function(object, newdata = NULL, t_grid = seq(0,1, 0.01), ...){
  # create x dataframe
  if(is.null(newdata)){
    x_data <- numeric()}
  else{
    x_data <- data.frame(newdata)
    x_data <- x_data[, all.vars(object$formula[[3]]), drop = FALSE]
  }
  x_data_old <- object$x_data[, all.vars(object$formula[[3]]), drop = FALSE]
  # create x model matrix
  if(ncol(x_data_old) == 0){
    x_design <- matrix(1)
  } else {
    factor_vars <- names(which(sapply(x_data_old, is.factor)))
    for(var in factor_vars){
      x_data[,var] <- factor(x_data[,var], levels = levels(x_data_old[,var[1]]))
    }
    x_design <- model.matrix(object$formula, cbind("data_curves" = 1, rbind(x_data, x_data_old)))
    if(!is.null(newdata)){
    x_design <- x_design[1:nrow(x_data),,drop = FALSE]}
  }

  lapply(1:nrow(x_design), function(k){
    x <- x_design[k,]
    pred_coefs <- Reduce("+", lapply(1:length(x), function(j) x[j]*object$coefs_list[[j]]))
    if(object$type == "smooth"){
      srv_curve <- function(t){
        t(make_reg_design(x_design_mat = NULL, t = t, type = "smooth",
                          knots = object$knots, closed = object$closed) %*% pred_coefs)
      }
      pred_data_curve <- as.data.frame(t(srvf_to_curve(t_grid, srv_curve)))
    } else {
      pred_data_curve <- get_points_from_srv(data.frame("t" = object$knots[-length(object$knots)],
                                                        pred_coefs))
    }
    if(object$closed){
      pred_data_curve <- project_curve_on_closed(pred_data_curve)
    }

    #compute translation
    offset <- apply(get_evals(cbind("t" = get_arc_length_param(pred_data_curve),
                                    pred_data_curve)), 2, mean)
    data.frame(t(t(pred_data_curve) - offset))
  })
}
