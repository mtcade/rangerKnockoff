local_grad.for_col <- function(
    j, # int, column
    forest,
    dframe, # Table of data
    bandwidths, # vector[ float ]; literal bandwidth for column `j`
    base.predictions, # vector[ float ] prediction from base forest without perturbing data
    verbose = 0
){ # -> vector[ float], length = dim(X)[1]
  if ( verbose > 0 & j %% 20 == 0 ){
    print(paste("# Taking local_grad.for_col", j, "of", dim(dframe)[2] ))
  }

  if ( length( bandwidths ) == 1 ){
    bandwidth <- bandwidths
  } else {
    bandwidth <- bandwidths[ j ]
  }
  if ( inherits( dframe[[j]], "factor") ){
    # Categorical
    X.j.levels <- levels( dframe[[j]] )
    levels.count <- length( X.j.levels )

    X.j.predictions <- matrix( data = 0, nrow = dim(dframe)[1], ncol = levels.count )
    for ( k in 1:levels.count ){
      X.test <- dframe
      X.test[[j]] <- X.j.levels[[k]]
      X.j.predictions[,k] <- predict( forest, X.test )$predictions
    }
    X.j.predictions.mins <- rep( 0, times = dim(dframe)[1] )
    X.j.predictions.maxs <- rep( 0, times = dim(dframe)[1] )
    for ( i in 1:dim(dframe)[1] ){
      X.j.predictions.mins[i] <- min( X.j.predictions[i,] )
      X.j.predictions.maxs[i] <- max( X.j.predictions[i,] )
    }
    local_grads <- X.j.predictions.maxs - X.j.predictions.mins
  } else {
    # Numeric
    X.test <- dframe

    # Choose randomly whether to go plus or minus bandwidth
    bandwidth.factor <- bandwidth * sample( c(-1,1), 1 )
    X.test[[j]] <- dframe[[j]] + bandwidth.factor

    local_grads <- predict( forest, data = X.test )$predictions - base.predictions
    local_grads <- local_grads / bandwidth.factor
  }

  return( local_grads )
}

stat.forest.local_grad.extras <- function(
    X,
    X_k,
    y,
    bandwidth = 1,
    bandwidth.exponent = 0.2,
    exponent = 2,
    in.parallel = FALSE,
    get.error = TRUE,
    get.oob_score = TRUE,
    verbose = 0,
    ...
){ # -> list("W","error","oob_score")
  # Row count must be equal
  stopifnot( all( dim(X) == dim(X_k) ) )

  # Get the matrix of gradients at each sample point
  # n factor for bandwidth; bandwidths[[ j ]] = sd(X[[j]])*bandwidth/n**n.factor
  n.factor <- dim(X)[1]**(bandwidth.exponent)

  vargs <- list( ... )
  if ( !( "respect.unordered.factors" %in% names(vargs) ) ){
    vargs[[ "respect.unordered.factors" ]] <- "partition"
  }

  X_all <- cbind( X, X_k )
  p_all <- dim(X_all)[2]

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
  base.predictions <- predict( forest, data = X_all )$predictions
  # Get the matrix of local grads

  # Calculate termwise bandwidths
  bandwidths <- rep( 0, times = p_all )
  for ( j in 1:p_all ){
    if ( inherits( X_all[[ j ]], "numeric" )  ){
      bandwidths[ j ] <- sd( X_all[[ j ]] )*bandwidth/n.factor
    }
  }

  if ( in.parallel ){
    library( parallel )
    no_cores <- detectCores()
    if ( verbose > 0 ){
      print(paste("# Running parallel with", no_cores, "cores"))
    }
    # Creating a cluster with the number of cores
    clust <- makeCluster(no_cores)

    local_grads.list <- parLapply(
      clust,
      1:p_all,
      fun = local_grad.for_col,
      forest = forest,
      dframe = X_all,
      bandwidths = bandwidths,
      base.predictions = base.predictions,
      verbose = verbose
    )
    print("parallel local_grads.list:")
    print( local_grads.list )
  } else {
    local_grads.list <- lapply(
      X = 1:p_all,
      FUN = local_grad.for_col,
      forest = forest,
      dframe = X_all,
      bandwidths = bandwidths,
      base.predictions = base.predictions,
      verbose = verbose
    )
  }

  local_grads.matrix <- do.call(
    cbind,
    local_grads.list
  )

  if ( in.parallel ){
    print("# local_grads.matrix:")
    print( local_grads.matrix )
  }

  # Transform to importances with abs and exponent
  importances <- colMeans( abs( local_grads.matrix )**exponent )

  # Calculate W from importances using differences, and pack results
  p_half <- p_all %/% 2
  results <- list(
    W = importances[1:p_half] - importances[(p_half+1):(p_all)]
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
    in.parallel = FALSE,
    verbose = 0,
    ...
){
  return(
    stat.forest.local_grad.extras(
      X = X,
      X_k = X_k,
      y = y,
      bandwidth = bandwidth,
      bandwidth.exponent = bandwidth.exponent,
      exponent = exponent,
      in.parallel = in.parallel,
      get.error = FALSE,
      get.oob_score = FALSE,
      verbose = verbose,
      ...
    )[[ "W" ]]
  )
}
