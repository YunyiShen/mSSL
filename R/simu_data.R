mprobit <- function(X,B,Sigma, Omega){
    unitdiag(Sigma,Omega)
    XB <- X %*% B
    E <- mvrnorm(nrow(X),rep(0,ncol(B)),Sigma)
    1.*(XB+E>=0)

}
