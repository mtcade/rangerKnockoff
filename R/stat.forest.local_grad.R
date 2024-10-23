#' Forest Local Gradient Importance
#'
#' @description Computes the variable importance from the mean absolute value of
#'  the local gradient approximation, and returns the difference between the original
#'  and the knockoff importances.
#' Handles categorical data by using the maximum prediction difference between categories.
#'
#' @param X n-by-p data.frame of original variables
#' @param X_k n-by-p data.frame of knockoff variables
#' @param y vector of length n, containing the response variables
#' @param bandwidth distance for estimating the numeric local gradients, multiplied by
#'  the standard deviation of the variable, divided by `n^bandwidth.power`
#' @param bandwidth.power exponent of sample size in the demoninator for bandwidth approximation
#'  of local gradient. Set to 0 to ignore sample size.
#' @param power Exponent for the absolute value of the estimated local gradient
#' @param ... Other values passed to [ranger::ranger]
#'
#' @return A vector of length p of W statistics
#'
#' @export

stat.forest.local_grad <- function(
  X,
  X_k,
  y,
  bandwidth = 1,
  bandwidth.power = 0.2,
  power = 2,
  ...
){
  stopifnot(
    all(
      dim(X_k) == dim(X)
    )
  )
  p <- dim(X)[2]
  # n factor for bandwidth; bandwidths[[ j ]] = std(X[[j]])*bandwidth/n**n.factor
  n.factor <- dim(X)[1]**(bandwidth.power)

  vargs <- list( ... )
  if ( !( "respect.unordered.factors" %in% names(vargs) ) ){
    vargs[[ "respect.unordered.factors" ]] <- "partition"
  }

  X_all <- cbind( X, X_k )

  # Calculate termwise bandwidths
  bandwidths <- rep( 0, times = 2*p )
  for ( j in 1:(2*p ) ){
    if ( inherits( X_all[[ j ]], "numeric" )  ){
      bandwidths[ j ] <- std( X_all[[ j ]] )*bandwidth/n.factor
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
    rbind,
    lapply(
      1:dim(X)[1],
      function ( i ) abs.local_grad(
        forest = forest,
        X = X_all[i,],
        bandwidths = bandwidths,
        power = power
      )
    )
  )

  importances <- colMeans( matrix.local_grad )

  return( importances[1:p] - importances[(p+1):2*p] )
}

abs.local_grad <- function(
  forest, # ranger::ranger
  X, # Vector of length 2p
  bandwidths, # vector[ float ]
  power # float
){
  stopifnot( length( bandwidths ) == length( X ) )
  X.abs.local_grad <- rep( 0, times = length( X ) )
  for ( j in 1:length(X) ){
    if (  inherits( X[[ j ]], "numeric") ){
      # Estimate local grad using bandwidths
      X.0 <- X
      X.0[[ j ]] <- X.0[[ j ]] - bandwidths[[ j ]]/2
      X.1[[ j ]] <- X
      X.1[[ j ]] <- X.1[[ j ]] + bandwidths[[ j ]]/2

      X.local_grad[[ j ]] <- ( predict( forest, X.1 ) - predict( forest, X.0 ) )/bandwidths[[ j ]]
    } else if ( inherits( X[[ j ]], "factor" ) ){
      print(paste("# Check X[[",j,"]]",sep = ""))
      print( levels( X[[j]] ) )
      stop("Check")
      predictions <- rep( 0, times = length( levels( X[[j]] ) ) )

      for ( k in 1:length( levels( X[[j]] ) ) ){
        X.k <- X
        X.k[[ j ]] <- levels( X[[j]] )[k]
        predictions[ k ] <- predict( forest, X.k )

        predictions.range <- range( predictions )
        X.abs.local_grad[[ j ]] <- predictions.range[2] - predictions.range[1]
      }

    } else {
      stop(paste("Unrecognized class X[[", j, "]] = ", class(X[[j]]), sep = "") )
    }

    # Take absolute value and power
    X.abs.local_grad[[ j ]] <- abs( X.abs.local_grad[[ j ]] )**power
  }

  return( X.abs.local_grad )
}
