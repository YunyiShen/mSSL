mprobit <- function(X,B,mu,Sigma, Omega){
    unitdiag(Sigma,Omega)
    XB <- X %*% B
    E <- mvrnorm(nrow(X),mu,Sigma)
    1.*(XB+E>=0)

}

cgprobit <- function(X,B,mu,Sigma, Omega){
  unitdiag(Omega,Sigma)
  XB <- X %*% B
  E <- mvrnorm(nrow(X),mu,Omega)
  res <- XB+E
  res <- res %*% Sigma
  1.*(res>=0)
}
