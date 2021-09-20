
isclose <- function(a,b, tol = 1e-10){
  abs(a-b)<=tol
}

#' Calculating several error measurement of Omega matrix
#' @name error_Omega
#' @description There are 11 measurement in total to be measured. 
#' @details \itemize{
#'  \item{TP}{:true positive rate}
#'    \item{TN}{:true negative rate}
#'    \item{FP}{:false positive rate}
#'    \item{FN}{:false negative rate}
#'    \item{SEN}{:sensitivity}
#'    \item{SPE}{:specificity}
#'    \item{PREC}{:precision}
#'    \item{ACC}{:accuracy }
#'    \item{F1}{:2 * TP/(2*TP + FP + FN)}
#'    \item{MCC}{:Matthew's Correlation Coefficie}
#'  \item{FROB:} {squared Frobenius error}
#' } 
#' @param Omega the estimated Omega
#' @param Omega_orig ground truth of Omega
#' @return a named vector of the 11 measurements

error_Omega <- function(Omega, Omega_orig)
{
  ut_Omega <- Omega[upper.tri(Omega, diag = FALSE)]
  ut_Omega_orig <- Omega_orig[upper.tri(Omega_orig, diag = FALSE)]
  
  perf <- rep(NA, times = 11)
  names(perf) <- c("TP", "TN", "FP", "FN", "SEN", "SPE", "PREC", "ACC",  "F1", "MCC", "FROB")
  tp <- sum(!isclose(ut_Omega, 0) & !isclose(ut_Omega_orig, 0))
  tn <- sum(isclose(ut_Omega, 0) & isclose(ut_Omega_orig , 0))
  fp <- sum(!isclose( ut_Omega, 0) & isclose(ut_Omega_orig , 0))
  fn <- sum(isclose(ut_Omega, 0) & !isclose(ut_Omega_orig, 0))
  
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
  perf["FROB"] <- sum( (Omega - Omega_orig)^2 )
  return(perf)
 } 
