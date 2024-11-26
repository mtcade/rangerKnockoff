abs.local_grad <- function(
  forest,
  X, # Table of data
  bandwidth, # float; literal bandwidth for column `j`
  exponent, # power of the local grads
  j # int, column
){ # -> vector[ float], length = dim(X)[1]

  if ( inherits( X[[j]], "factor") ){
    # Categorical
    X.j.levels <- levels( X[[j]] )
    levels.count <- length( X.j.levels )

    X.j.predictions <- matrix( data = 0, nrow = dim(X)[1], ncol = levels.count )
    for ( k in 1:levels.count ){
      X.test <- X
      X.test[[j]] <- X.j.levels[[k]]
      X.j.predictions[,k] <- predict( forest, X.test )$predictions
    }
    X.j.predictions.mins <- rep( 0, times = dim(X)[1] )
    X.j.predictions.maxs <- rep( 0, times = dim(X)[1] )
    for ( i in 1:dim(X)[1] ){
      X.j.predictions.mins[i] <- min( X.j.predictions[i,] )
      X.j.predictions.maxs[i] <- max( X.j.predictions[i,] )
    }
    local_grads <- X.j.predictions.maxs - X.j.predictions.mins
  } else {
    # Numeric
    X.test <- X
    X.test[[j]] <- X[[j]] + bandwidth*0.5

    local_grads <- predict( forest, data = X.test )$predictions

    X.test[[j]] <- X[[j]] - bandwidth*0.5

    local_grads <- local_grads - predict( forest, data = X.test )$predictions

    local_grads <- abs(local_grads) / bandwidth
  }

  return( local_grads**exponent )
}

stat.forest.local_grad.extras <- function(
  X,
  X_k,
  y,
  bandwidth = 1,
  bandwidth.exponent = 0.2,
  exponent = 2,
  get.error = TRUE,
  get.oob_score = TRUE,
  ...
){ # -> list("W","error","oob_score")
  stopifnot(
    all(
      dim(X_k) == dim(X)
    )
  )
  p <- dim(X)[2]
  # n factor for bandwidth; bandwidths[[ j ]] = sd(X[[j]])*bandwidth/n**n.factor
  n.factor <- dim(X)[1]**(bandwidth.exponent)

  vargs <- list( ... )
  if ( !( "respect.unordered.factors" %in% names(vargs) ) ){
    vargs[[ "respect.unordered.factors" ]] <- "partition"
  }
  vargs[[ "oob.error" ]] <- get.oob_score

  X_all <- cbind( X, X_k )

  # Calculate termwise bandwidths
  bandwidths <- rep( 0, times = 2*p )
  for ( j in 1:(2*p ) ){
    if ( inherits( X_all[[ j ]], "numeric" )  ){
      bandwidths[ j ] <- sd( X_all[[ j ]] )*bandwidth/n.factor
    }
  }

  forest <- do.call(
    ranger::ranger,
    c(
      list(
        x = X_all,
        y = y
      ),
      vargs
    )
  )
  matrix.local_grad <- do.call(
    cbind,
    lapply(
      1:(2*p),
      function( x ) abs.local_grad(
        forest = forest,
        X = X_all,
        bandwidth = bandwidths[ x ],
        exponent = exponent,
        j = x
      )
    )
  )

  importances <- colMeans( matrix.local_grad )

  # Pack results
  results <- list(
    W = importances[1:p] - importances[(p+1):(2*p)]
  )
  if ( get.error ){
    results[[ "error" ]] <- var( predict( forest, X_all)$predictions - y )
  }
  if ( get.oob_score ){
    results[[ "oob_score" ]] <- forest[[ "prediction.error" ]]
  }

  return( results )
}

#' Forest Local Gradient Importance
#'
#' @description Computes the variable importance from the mean absolute value of
#'  the local gradient approximation, and returns the difference between the original
#'  and the knockoff importances.
#' Handles categorical data by using the maximized prediction difference between categories.
#'
#' @param X n-by-p data.frame of original variables
#' @param X_k n-by-p data.frame of knockoff variables
#' @param y vector of length n, containing the response variables
#' @param bandwidth distance for estimating the numeric local gradients, multiplied by
#'  the standard deviation of the variable, divided by `n^bandwidth.exponent`
#' @param bandwidth.exponent exponent of sample size in the demoninator for bandwidth approximation
#'  of local gradient. Set to 0 to ignore sample size.
#' @param exponent Exponent for the absolute value of the estimated local gradient
#' @param ... Other values passed to [ranger::ranger]
#'
#' @return A vector of length p of W statistics
#'
#' @examples
#' stat.forest.local_grad( X, X_k, y, bandwidth = 0.5, exponent = 1 )
#'
#'
#' @export

stat.forest.local_grad <- function(
  X,
  X_k,
  y,
  bandwidth = 1,
  bandwidth.exponent = 0.2,
  exponent = 2,
  ...
){
  return(
    stat.forest.local_grad.extras(
      X = Xk,
      X_k = X_k,
      y = y,
      bandwidth = bandwidth,
      bandwidth.exponent = bandwidth.exponent,
      exponent = exponent,
      get.error = FALSE,
      get.oob_score = FALSE,
      ...
    )[[ "W" ]]
  )
}
