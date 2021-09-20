# Stability plot


# This only makes sense for objects resulting from mSSL_dpe
# Look at the output for early_term.
# Early_term == 1 means that we S was ill-conditioned
# current_index = s + t*L (in C++ with 0-index). In R it would be (s-1) + (t-1)*L + 1

# Loop over the results in B_path and Omega_path
# Identify the points that have same support as the terminal point. If it does, color the point blue, give it pch = 3 (+). If it does not, color the point black, give it pch = 16 (filled circle). If it is unstable, color the point red, give it pch = 4 (x)\n


# In the source code for mSSL_dpe we loop over the lambdas and xi's
# for(s in 1:L){
  # for(t in 1:L){
    # lambda0 = lambda_0[s]; xi0 = xi_s[t];
#  }
#}

#' Stability plot of a dpe result
#' @description This only makes sense for objects resulting from dynamic posterior exploration. Identify the points that have same support as the terminal point.
#' @param fit_cgSSL_dpe output object from cgSSL that uses condexp=FALSE
#' @details 
#' Loop over the results in B_path and Omega_path
#' Identify the points that have same support as the terminal point. If it does, color the point blue, give it pch = 3 (+). If it does not, color the point black, give it pch = 16 (filled circle). If it is unstable, color the point red, give it pch = 4 (x)

plot_support <- function(fit_cgSSL_dpe){
  
  lambda0 <- fit_cgSSL_dpe$lambda0
  xi0 <- fit_cgSSL_dpe$xi0
  L <- length(lambda0)
  supp_B_final <- which(fit_cgSSL_dpe$B != 0)
  supp_Omega_final <- which(fit_cgSSL_dpe$Omega != 0)
  
  par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0), cex.axis = 0.8, cex.lab = 0.8, cex.axis = 0.8)
  plot(1, type = "n", xlim = range(lambda0), ylim = range(xi0), xlab = expression(lambda[0]), ylab = expression(xi[0]), main = "")
  for(s in 1:L){
    for(t in 1:L){
      index <- (s - 1) + L*(t-1) + 1
      supp_B <- which(fit_cgSSL_dpe$B_path[,,index] != 0)
      supp_Omega <- which(fit_cgSSL_dpe$Omega_path[,,index] != 0)
      
      if(fit_cgSSL_dpe$early_term[index] == 1){
        points(lambda0[s], xi0[t], col = 'red', pch = 4)
      }
      if(identical(supp_B, supp_B_final) & identical(supp_Omega, supp_Omega_final) & fit_cgSSL_dpe$early_term[index] == 0){
        points(lambda0[s],xi0[t], col = 'blue', pch = 3)
      } else if(fit_cgSSL_dpe$early_term[index] == 0){
        points(lambda0[s],xi0[t], col = 'black', pch = 16)
      } else if(fit_cgSSL_dpe$early_term[index] == 1){
        points(lambda0[s],xi0[t], col = 'red', pch = 4)
      }
    }
  }
  
  # Add legend at the top or on the right
  # Red x: numerically unstable S
  # Blue plus: stable support.
  
  
}
