library( rangerKnockoff )
data_X <- listRDA("/Users/evanmason/Documents/UCR/r_custom/rangerKnockoff/data/X.rda")
X <- data_X[["X"]]

beta <- 0.5
y <- X[["0"]] - 3.76*(1/(X[["5"]]^2 +1))
y <- y + 1.42*cos(2*pi*X[["28"]]) + 2.86*sqrt(abs(X[["84"]]))
y <- y - 0.7*(X[["95"]]^2)
y <- y + 1*(X[["96"]] == "1") - 1*(X[["110"]] == "1")
y <- y - 1 *(X[["127"]] == "1") + 1*(X[["127"]] == "2")
y <- beta*y + rnorm( n=dim(X)[1] )

Xk <- rangerKnockoff::create.forest.conditional_residuals(
  X = X
)

W <- rangerKnockoff::stat.forest.local_grad( X = X, X_k = Xk, y = y )

knockoff_threshold <- knockoff::knockoff.threshold( W, fdr = 0.2 )

selected <- W >= knockoff_threshold

relevant <- c(1,6,11,16,20,21,25,33)

if ( sum( selected ) <= 0  ){
  power <- 0
  fdr <- 0
} else {
  power <- sum( selected[relevant] )/length( relevant )
  fdr <- 1 - sum( selected[relevant] )/sum( selected )
}

cat("Power:", power,"\n")
cat("fdr:", fdr,"\n")