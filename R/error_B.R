isclose <- function(a,b, tol = 1e-15){
  abs(a-b)<=tol
}

#' @title Calculating several error measurement of B matrix
#' @name error_B
#' @description There are 11 measurement in total to be measured.
#' @param B the estimated B
#' @param B_orig ground truth of B
#' @return a named vector of the 11 measurements
#' @details 
#' \itemize{
#'    \item{TP}{:true positive rate}
#'    \item{TN}{:true negative rate}
#'    \item{FP}{:false positive rate}
#'    \item{FN}{:false negative rate}
#'    \item{SEN}{:sensitivity}
#'    \item{SPE}{:specificity}
#'    \item{PREC}{:precision}
#'    \item{ACC}{:accuracy }
#'    \item{F1}{:2 * TP/(2*TP + FP + FN)}
#'    \item{MCC}{:Matthew's Correlation Coefficie}
#'    \item{MSE}{:mean squared error}
#' }
#' 
error_B <- function(B, B_orig){
  
    perf <- rep(NA, times = 11)
    names(perf) <- c("TP", "TN", "FP", "FN", "SEN", "SPE", "PREC", "ACC",  "F1", "MCC", "MSE")
    tp <- sum(!isclose(B, 0) & !isclose(B_orig, 0))
    tn <- sum(isclose(B, 0) & isclose(B_orig , 0))
    fp <- sum(!isclose( B, 0) & isclose(B_orig , 0))
    fn <- sum(isclose(B, 0) & !isclose(B_orig, 0))
    
    
    sen <- tp/(tp + fn)
    spe <- tn/(tn + fp)
    prec <- tp/(tp + fp)
    acc <- (tp + tn)/(tp + tn + fp + fn)
    f1 <- 2 * tp/(2*tp + fp + fn)
    
    mcc.n <- tp + tn + fp + fn
    mcc.s <- (tp + fn)/mcc.n
    mcc.p <- (tp + fp)/mcc.n
    mcc <- (tp/mcc.n - mcc.s*mcc.p)/sqrt(mcc.p * mcc.s * (1 - mcc.p) * (1 - mcc.s))
    perf["TP"] <- tp 
    perf["TN"] <- tn 
    perf["FP"] <- fp 
    perf["FN"] <- fn 
    perf["SEN"] <- sen
    perf["SPE"] <- spe
    perf["PREC"] <- prec
    perf["ACC"] <- acc
    perf["F1"] <- f1
    perf["MCC"] <- mcc
    perf["MSE"] <- mean((B - B_orig)^2)
  return(perf)

}

log_l2 <- function(a,b){
    log(sum((a-b)^2))
}