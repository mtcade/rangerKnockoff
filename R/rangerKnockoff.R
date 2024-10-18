#' Knockoffs with conditional residuals
#'
#' @param X data.frame requiring knockoffs
#' @param method `c("second_order","fixed","none")`
#' Procedure for creating numeric conditional residuals.
#' Ignored if `residuals_function` is provided.
#'  "second_order": Gaussian second order knockoffs; see `knockoff::create.second_order`
#'  "fixed": Fixed knockoffs; see `knockoff:create.fixed`
#'  "none": Only calculate conditional expectations. Useful for creating knockoffs with other software.
#' @param residuals_method `c("asdp","equi","sdp")`
#' Procedure for `method` of `"second_order"` or `"fixed"` when calculating covariance.
#' Default is "asdp" for `method = "second_order"`, "sdp" for `method = "fixed"`
#'  "asdp": See `knockoff::create_asdp`
#'  "equi": See `knockoff::create_equi`
#'  "sdp": See `knockoff::create_sdp`
#' @param residuals_function If provided, uses any custom function to create numeric residuals.
#'  Must take an arbitrary sized numeric data.frame as an input, and return a numeric data.frame of the same size.
#' @param sigma Optional covariance when `method = "fixed"`
#'  See `knockoff::create_fixed`
#' @param randomize Optional bool when `method = "fixed"`
#'  See `knockoff::create_fixed`
#' @param shrink Optional bool when `method = "fixed"`
#'  See `knockoff::create_fixed`
#' @param ... Other values passed to `ranger::ranger`
#' @return A knockoff data.frame with the same size, column names, and column types as `X`
#' @export
create.forest.conditionalResiduals <- function(
  X, # n-by-p matrix of original variables
  method = "second_order", # How to create numeric conditional residual knockoffs
  # "second_order", "fixed": from `knockoff` package
  # "none": Only conditional expectations for numeric, still pick categories
  #   Used mostly for exporting to python for more sophisticated methods
  residuals_method = NULL, # c("asdp","equi","sdp")
  residuals_function = NULL, # Custom function to turn a dattable of residuals into knockoff residuals
  #   `method` and `residuals_method` are ignored if this is non null
  sigma = NULL, # For method = "fixed"
  randomize = F, # For method = "fixed"
  shrink = F, # For method = "second_order"
  ... # inputs for ranger
){
  # -- Wrangle inputs
  vargs <- list( ... )
  if ( !( "respect.unordered.factors" %in% names(vargs) ) ){
    vargs[[ "respect.unordered.factors" ]] <- "partition"
  }

  #
  p = ncol(X) # int
  X.names <- names(X)

  forestPredictions <- lapply(
    1:p,
    FUN = function( j ) do.call(
      get.forestPredictions,
      c(
        list(
          X = X,
          dependent.variable.name = X.names[[ j ]]
        ),
        vargs
      )
    )
  ) # list[ vector[numeric] | matrix[numeric] ]

  X.classes <- lapply(
    X,
    class
  ) # list[ "numeric" | "factor" ]

  # Categorical Variables
  forestRandomCategories <- lapply(
    X = forestPredictions[
      X.classes == "factor"
    ],
    FUN = function(x) categoryChooser( x )
  ) # list[ vector[integer] ]

  # Make the numeric data predictions into a matrix,
  # To calculate the conditional residuals
  conditionalExpectationsNumeric <- do.call(
    cbind,
    forestPredictions[
      X.classes == "numeric"
    ]
  ) # matrix[ numeric ]

  if ( !is.null(residuals_function) ){
    conditionalResidualsNumeric <- matrix(
      data = 0.0,
      nrow = dim( conditionalExpectationsNumeric )[1],
      ncol = dim( conditionalExpectationsNumeric )[2]
    )

    .numericIterator <- 1
    for ( j in 1:p ){
      if ( X.classes[[j]] == "numeric" ){
        conditionalResidualsNumeric[
          ,.numericIterator
        ] <- X[,j] - conditionalExpectationsNumeric[
          ,.numericIterator
        ]
        .numericIterator <- .numericIterator + 1
      }
    }

    conditionalResidualKnockoffs <- residuals_function( conditionalResidualsNumeric )
  }
  else if ( method == "none" ){
    # Use only conditional expectations
    .categoryIterator <- 1
    .numericIterator <- 1

    Xk <- X
    for ( j in 1:p ){
      if ( X.classes[[j]] == "factor" ){
        # Convert chosen indices into the levels from X for categories
        Xk[[j]] <- levels( X[[j]] )[ forestRandomCategories[[ .categoryIterator ]] ]
        .categoryIterator <- .categoryIterator + 1
      } else {
        # Use conditional residual knockoffs for numerics
        Xk[[j]] <- conditionalExpectationsNumeric[,.numericIterator]
        .numericIterator <- .numericIterator + 1
      }
    }

    stopifnot( .categoryIterator + .numericIterator - 2 == p )
  }
  else {
    # Calculate conditional residuals and their knockoffs
    conditionalResidualsNumeric <- matrix(
      data = 0.0,
      nrow = dim( conditionalExpectationsNumeric )[1],
      ncol = dim( conditionalExpectationsNumeric )[2]
    )

    .numericIterator <- 1
    for ( j in 1:p ){
      if ( X.classes[[j]] == "numeric" ){
        conditionalResidualsNumeric[
          ,.numericIterator
        ] <- X[,j] - conditionalExpectationsNumeric[
          ,.numericIterator
        ]
        .numericIterator <- .numericIterator + 1
      }
    }

    # make the knockoffs
    if ( method == "fixed" ){
      if ( is.null( residuals_method ) ){
        residuals_method <- "sdp"
      }
      residuals_method <- match.arg( residuals_method, c("sdp","equi") )

      conditionalResidualKnockoffs <- knockoff::create.fixed(
        X = conditionalResidualsNumeric,
        method = residuals_method,
        sigma = sigma,
        randomize = randomize
      )
    } else if ( method == "second_order" ){
      if ( is.null( residuals_method ) ){
        residuals_method <- "asdp"
      }
      residuals_method <- match.arg( residuals_method, c("asdp","equi","sdp") )

      conditionalResidualKnockoffs <- knockoff::create.second_order(
        X = conditionalResidualsNumeric,
        method = residuals_method,
        shrink = shrink
      )
    } else {
      stop( paste( "Unrecognized method:", method ) )
    }

    # Make the knockoff data by taking X and replacing columns from knockoffs
    .numericIterator <- 1
    .categoryIterator <- 1

    Xk <- X
    for ( j in 1:p ){
      if ( X.classes[[j]] == "factor" ){
        # Convert chosen indices into the levels from X for categories
        Xk[[j]] <- levels( X[[j]] )[ forestRandomCategories[[ .categoryIterator ]] ]
        .categoryIterator <- .categoryIterator + 1
      } else {
        # Use conditional residual knockoffs for numerics
        Xk[[j]] <- conditionalExpectationsNumeric[,.numericIterator] + conditionalResidualKnockoffs[,.numericIterator]
        .numericIterator <- .numericIterator + 1
      }
    }

    stopifnot( .categoryIterator + .numericIterator - 2 == p )
  }

  colnames(Xk) <- paste( colnames(X), "~", sep = "" )

  return( Xk )
}

#' Knockoffs with SCIP
#'
#' @export
create.forest.SCIP <- function(
  X, # n-by-p matrix of original variables
  method, # "normal","permute",
  ... # inputs for ranger
){
  # -- Wrangle inputs

  method <- match.arg( method, c("normal","permute") )
  # Arguments; set defaults different from ranger
  vargs <- list( ... )
  if ( !( "respect.unordered.factors" %in% names(vargs) ) ){
    vargs[[ "respect.unordered.factors" ]] <- "partition"
  }

  #
  p = ncol(X) # int
  X.names <- names(X)

  Xk <- do.call(
    get.forestPredictions.forColumn,
    c(
      list( X = X, column = 1, method = method, Xk = NULL),
      vargs
    )
  )

  stopifnot( dim(Xk) > 1 )

  for ( j in 2:dim(X)[2] ){
    Xk.j <- do.call(
      get.forestPredictions.forColumn,
      c(
        list( X = X, column = j, method = method, Xk = Xk),
        vargs
      )
    )

    Xk <- cbind( Xk, Xk.j )
  }


  colnames( Xk ) <- paste( colnames( X ), "~", sep = "" )

  return( Xk )
}

# -- Helper Functions

get.forestPredictions <- function(
    X, # n-by-p
    dependent.variable.name, # int
    ... # ranger options; do not use 'probability' or 'dependent.variable.name'
){
  X.class <- class( X[[ dependent.variable.name ]] )

  if ( X.class == "numeric" ){
    probability <- F
  } else if ( X.class == "factor" ){
    probability <- T
  } else {
    stop(paste("Unrecognized class( X[[ dependent.variable.name ]] ):", X.class))
  }

  forest <- ranger::ranger(
    data = X,
    dependent.variable.name = dependent.variable.name,
    probability = probability,
    ...
  )

  predictions <- forest$predictions
  return( predictions )
}

get.forestPredictions.forColumn <- function(
  X, # n by p data frame
  column, # int; index
  method = "normal", # "normal", "permute"
  Xk = NULL,
  ...
){
  # Predict Xk[[ column ]] based on X[1:n,-column] and Xk[ 1:n, 1:(column-1) ]
  # Generate new Xk based on predictions and residuals_method for numeric X[[ column ]],
  #    pick based on probabilities for factor X[[ column ]]
  n <- dim( X )[1]
  stopifnot( dim(Xk)[1] == n )

  method <- match.arg( method, c("normal","permute") )

  if ( column == 1 ){
    # No Xk at all
    X.all <- X
  }
  else if ( column == 2 ){
    # Xk will be a vector, not a dataframe
    X.all <- cbind( X, Xk )
  }
  else {
    # Xk is a data frame
    Xk.temp <- Xk[ 1:n, 1:(column-1) ]
    X.all <- cbind( X, Xk.temp )
  }

  if ( class(X[[ column ]]) == "numeric" ){
    forest <- ranger::ranger(
      x = X.all[1:n, -column],
      y = X.all[1:n, column],
      probability = F,
      ...
    )

    residuals <- X[[ column ]] - forest$predictions

    if ( method == "normal" ){
      residuals.Xk.j <- rnorm(
        n = n,
        mean = 0,
        sd = sd( residuals )
      )
    } else if ( method == "permute" ){
      residuals.Xk.j <- sample( residuals )
    } else {
      stop( paste( "Unrecognized residuals_method:", residuals_method ) )
    }

    Xk.j <- forest$predictions + residuals.Xk.j
  }
  else if ( class(X[[ column ]]) == "factor" ){
    forest <- ranger::ranger(
      x = X.all[1:n, -column],
      y = X.all[1:n, column],
      probability = T,
      ...
    )

    Xk.j <- as.factor(
      categoryChooser( forest$predictions )
    )
    levels( Xk.j ) <- levels( X[[ column ]] )
  }
  else {
    stop(
      paste(
        "Unrecognized class( X[[ column ]] ):", class( X[[ column ]] )
      )
    )
  } #/switch class( X[[ column ]] )

  return( Xk.j )
}

categoryChooser <- function( X ){
  # For a matrix X, choose from the number of columns with the rowwise
  #   probabilities. One chosen for each row as a column integer.
  # HINT: use the column integer to choose the appropriate factor level
  k <- dim( X )[2]
  choices <- rep( 0, times = dim(X)[1] )

  choices <- unlist(
    lapply(
      1:dim(X)[1],
      function(x) sample( k, size = 1, prob = X[x,] )
    )
  )

  return( choices )
}
