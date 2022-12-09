#' Plot method for planar elastic regression models
#' @description Plots objects of class \code{elastic_reg_model}.
#' @param x object of class \code{elastic_reg_model},
#' usually a result of a call to \code{\link{fit_elastic_regression}}
#' @param asp numeric, giving the aspect ratio of the two coordinates,
#' see \code{\link{plot.window}} for details.
#' @param col color of the predicted curves.
#' @param ... further plotting parameters.
#' @importFrom graphics plot lines legend
#' @importFrom grDevices colorRampPalette
#' @return No value
#' @export
#'
#' @seealso For examples see documentation of \code{\link{fit_elastic_regression}}.

plot.elastic_reg_model <- function(x, asp = 1, col = "red", ...){
  if(!(ncol(x$coefs_list$beta_0) == 2)){
    stop("Plotting option only for planar curves!")
  }
  col_names <- colnames(x$coefs_list$beta_0)
  data_curves <- lapply(x$data_curves, function(data_curve){
    offset <- apply(get_evals(cbind("t" = get_arc_length_param(data_curve),
                                    data_curve)), 2, mean)
    data.frame(t(t(data_curve) - offset))
  })
  data_curves <- lapply(data_curves, function(data) data[,col_names])
  data_all <- do.call("rbind", data_curves)

  #find covariates for predictions
  x_data <- x$x_data[, names(x$x_data) %in% all.vars(x$formula[[3]]), drop = FALSE]
  model_matrix <- model.matrix(x$formula, cbind("data_curves" = 1, x_data))
  factor_vars <- names(attributes(model_matrix)$contrasts)
  dummydata <- x_data[1,, drop = FALSE]
  dummydata[, !(colnames(dummydata) %in% factor_vars)] <-
    round(apply(x_data[, !(colnames(dummydata) %in% factor_vars), drop = FALSE], 2, mean))
  dummydata[, colnames(dummydata) %in% factor_vars] <-
    apply(x_data[, factor_vars, drop = FALSE], 2, function(x) names(which.max(table(x))))

  for(i in 1:ncol(x_data)){
    #empty plot
    plot(NULL, xlim = range(data_all[,1]), ylim = range(data_all[,2]), xlab = col_names[1],
         ylab = col_names[2], asp = asp)
    #plot data
    invisible(lapply(data_curves, lines, col = "gray"))

    if(ncol(x_data) == 0){
      mean <- predict(x)[[1]]
      lines(mean, col = col, lwd = 2)
    } else if(colnames(x_data)[i] %in% factor_vars){
      newdata <- dummydata[rep(1,2),, drop = FALSE]
      factor_levels <- levels(factor(x_data[,i]))
      newdata[,i] <- factor_levels
      predictions <- predict(x, newdata = newdata)

      pred_colours <- colorRampPalette(c(col, "black"))(length(predictions))

      invisible(lapply(1:length(predictions), function(i){
        lines(predictions[[i]], col = pred_colours[[i]])
      }))
      legend("topright", col = pred_colours, legend = factor_levels, lty = 1)
    } else {
      newdata <- dummydata[rep(1,11),, drop = FALSE]
      range_x <- range(x_data[,i])
      mean_x <- mean(x_data[,i])
      delta <- min(mean_x - range_x[1], range_x[2] - mean_x)
      newdata[,i] <- mean_x + delta*(-5:5)/5
      predictions <- predict(x, newdata = newdata)

      pred_colours <- colorRampPalette(c("white", col, "black"))(length(predictions)*1.5)
      #plot predictions
      invisible(lapply(1:length(predictions), function(i){
        lines(predictions[[i]], col = pred_colours[[i + 3]])
      }))
      #emphasis mean covariate prediction
      lines(predictions[[6]], col = col, lwd = 2)
      legend("topright", col = pred_colours[c(1,6,11) + 3],
             legend = signif(newdata[c(1,6,11), i], 2), lty = 1,
             lwd = c(1,2,1))
    }
  }
}

