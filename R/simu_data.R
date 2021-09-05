mprobit <- function(X,B,mu,Sigma, Omega){
    unitdiag(Sigma,Omega)
    XB <- X %*% B
    E <- mvrnorm(nrow(X),mu,Sigma)
    1.*(XB+E>=0)

}
